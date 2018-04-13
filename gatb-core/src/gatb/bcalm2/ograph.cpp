#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>


/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */

// The code used to be optimized for uint128, but then i switched everything to GATB's LargeInt. It's likely unoptimized now.

using namespace std;

namespace gatb { namespace core { namespace debruijn { namespace impl  {

static string reverseinplace(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
		str[i] ^= 4;
		str[j] ^= 4;
		if ((str[i]&3) != 3){str[i]^= 17;}
		if ((str[j]&3) != 3){str[j]^= 17;}
		swap(str[i],str[j]);
	}
	if(str.size()%2==1){
		str[j] ^= 4;
		if ((str[j]&3) != 3){str[j]^= 17;}
	}
	return str;
}


static void reverseinplace2(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
		str[i] ^= 4;
		str[j] ^= 4;
		if ((str[i]&3) != 3){str[i]^= 17;}
		if ((str[j]&3) != 3){str[j]^= 17;}
		swap(str[i],str[j]);
	}
	if(str.size()%2==1){
		str[j] ^= 4;
		if ((str[j]&3) != 3){str[j]^= 17;}
	}
}


static uint chartoint(char c){
	char d = (c >> 1) & 3;
	if (d > 1)
		d ^= 1;
	return d;
}


static inline bool isNumber(char c){return (c<64);}


template<size_t span>
typename graph3<span>::kmerType graph3<span>::beg2int128(const string& str){
	typename graph3<span>::kmerType resBeg;
    resBeg.setVal(0);
	for(uint i(0);i<k;++i){
		resBeg= resBeg << 2;
		resBeg= resBeg + chartoint(str[i]);
	}
	return resBeg;
}


template<size_t span>
typename graph3<span>::kmerType graph3<span>::beg2int128rc(const string& str){
	typename graph3<span>::kmerType res;
    res.setVal(0);
	for(int i(k-1);i>=0;i--){
		res=res<<2;
		res=res+3-chartoint(str[i]);
	}
	return res;
}


template<size_t span>
typename graph3<span>::kmerType graph3<span>::end2int128rc(const string& str){
	typename graph3<span>::kmerType res;
    res.setVal(0);
	for(int i(k-1);i>=0;i--){
		res = res << 2;
		res = res + (3-chartoint(str[str.size()-k+i]));
	}
	return res;
}


template<size_t span>
typename graph3<span>::kmerType graph3<span>::end2int128(const string& str){
	typename graph3<span>::kmerType resEnd;
    resEnd.setVal(0);
	for(uint i(0);i<k;++i){
		resEnd= resEnd << 2;
		resEnd= resEnd + chartoint(str[str.size()-k+i]);
	}
	return resEnd;
}


template<size_t span>
// is that revcomp? if so, should use gatb's
typename graph3<span>::kmerType graph3<span>::rcb(typename graph3<span>::kmerType min){
	typename graph3<span>::kmerType resrcb;
	typename graph3<span>::kmerType three;
    three.setVal(3);
    resrcb.setVal(0);
	for(uint i(0); i<k;++i){
		resrcb= resrcb + ((three-(min%4))<<(2*(k-1-i)));
		min= min >> 2;
	}
	return resrcb;
}


template<size_t span>
void graph3<span>::compaction(uint iL,  uint iR,typename graph3<span>::kmerType kmmer){
	if(iR!=iL){
		typename graph3<span>::kmerType RC=rcb(kmmer);
		//uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
		bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
		if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]),kmmer);}
		if(b1){return compaction(stoi(unitigs[iL]),iR,kmmer);}
		if(b2){return compaction(iL,stoi(unitigs[iR]),kmmer);}
		//~ cout<<unitigs[iR]<<"\n";
		//~ cout<<unitigs[iL]<<"\n";
		typename graph3<span>::kmerType beg1;//(beg2int128(unitigs[iL])); // that kind of initialization isn't supported in LargeInt.
        beg1.setVal(beg2int128(unitigs[iL]));
		typename graph3<span>::kmerType end2;//(end2int128(unitigs[iR]));
        end2.setVal(end2int128(unitigs[iR]));

		if(beg1==end2 and (end2==kmmer or end2==RC)){
		//~ if(beg1==end2 ){
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
            indexed_right[iR]=indexed_right[iL];
            connected_right[iR]=connected_right[iL];
            compact_abundances(iR,iL);
			return;
		}

		typename graph3<span>::kmerType endrc2;//(beg2int128rc(unitigs[iR]));
        endrc2.setVal(beg2int128rc(unitigs[iR]));
		if(beg1==endrc2 and (beg1==kmmer or beg1==RC)){
		//~ if(beg1==endrc2 ){
			reverseinplace2(unitigs[iR]);
            indexed_left[iR]=indexed_right[iR];
            connected_left[iR]=connected_right[iR];

			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
            indexed_right[iR]=indexed_right[iL];
            connected_right[iR]=connected_right[iL];
            compact_abundances(iR,iL,true,false);
			return;
		}

		typename graph3<span>::kmerType beg2;//(rcb(endrc2));
        beg2.setVal(rcb(endrc2));
		typename graph3<span>::kmerType end1;//(end2int128(unitigs[iL]));
        end1.setVal(end2int128(unitigs[iL]));
		if(end1==beg2 and (end1==kmmer or end1==RC)){
		//~ if(end1==beg2 ){
			unitigs[iL]+=(unitigs[iR].substr(k));
			unitigs[iR]=to_string(iL);

            indexed_right[iL]=indexed_right[iR];
            connected_right[iL]=connected_right[iR];
            compact_abundances(iL,iR);
			return;
		}

		typename graph3<span>::kmerType begrc2;//(rcb(end2));
        begrc2.setVal(rcb(end2));
/*            std::cout << "a : " << rcb(end2).toString(31) << std::endl;
            std::cout << "a=b: " << begrc2.toString(31) << std::endl;*/ // manifestation of a bug when LargeInt constructors are removed

		if(end1==begrc2 and (end1==kmmer or end1==RC)){
		//~ if(end1==begrc2 ){
			unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
			unitigs[iR]=to_string(iL);

            indexed_right[iL]=indexed_left[iR];
            connected_right[iL]=connected_left[iR];
            compact_abundances(iL,iR, false, true);
			return;
		}
	}
}


template<size_t span>
void graph3<span>::compact_abundances(uint i1, uint i2, bool reverse_first, bool reverse_second){
    // aliases just to be able to copypaste stackoverflow code!
    std::vector<unsigned int>& v = unitigs_abundances[i1];
    std::vector<unsigned int>& v_prime = unitigs_abundances[i2];
    if (reverse_first)
        std::reverse(v.begin(),v.end());
    v.reserve(v.size() + distance(v_prime.begin(),v_prime.end()));
    if (reverse_second)
        std::reverse(v_prime.begin(),v_prime.end());
    v.insert(v.end(),v_prime.begin(),v_prime.end());
    v_prime.clear();
}

template<size_t span>
inline void graph3<span>::update_connected(kmerIndiceT<span> &ki)
{
    if (ki.position == SEQ_LEFT)
        connected_left[ki.indice] = true;
    if (ki.position == SEQ_RIGHT)
        connected_right[ki.indice] = true;
}

/*
 * this function is the core one that decides what compactions need to be made
 */
template<size_t span>
void graph3<span>::debruijn(){
	sort(left.begin(),left.end(),comparator<span>());
	sort(right.begin(),right.end(),comparator<span>());
	uint64_t iL(0),iR(0),sizeLeft(left.size()),sizeRight(right.size());
    typename graph3<span>::kmerType minusone;
    minusone.setVal(-1);
	left.push_back({0,minusone, SEQ_LEFT}); // dummy kmer so that we dont need to check bounds.. clever..
	right.push_back({0,minusone, SEQ_LEFT});
    uint debug_index = 0;

    for (uint32_t i = 0; i< indiceUnitigs; i++)
    {
        connected_left[i]  = false;
        connected_right[i] = false;
    }

    kmerIndiceT<span> kL,kR;
    std::vector<std::pair<uint,uint>> to_compact;
    // in this pass we just flag the pairs to compact.
    // before, we used to compact on the fly, but now i want to have proper connection info for all unitigs prior to compaction
	while(iL < sizeLeft && iR < sizeRight){
		kL=left[iL];
		kR=right[iR];

        //~ std::cout << " kl / kR " << kL.indice << " " << kR.indice << " " << kL.kmmer << " " << kR.kmmer << " unitigs " << unitigs[kL.indice] << " " << unitigs[kR.indice] << std::endl;

		if(kL.kmmer==kR.kmmer){
            //~ if (debug_index > 0) if (kL.indice == debug_index || kR.indice == debug_index ) std::cout << " identical, kl / kR " << kL.indice << " " << kR.indice << " unitigs " << unitigs[kL.indice] << " " << unitigs[kR.indice] << " positions "  << kL.position << " " << kR.position << std::endl;
            //~ if(isNumber (unitigs[kL.indice][0])){
			//~ }
            if(not kL.indice==kR.indice){
				update_connected(kL);
				update_connected(kR);
			}

            // found the same (k-1)-mer in the left and right array, it means that two sequences end with those and could be potentially compacted
			bool go(true);
			++iL;++iR;
			if(left[iL].kmmer==kL.kmmer){
				go=false;
				if(not left[iL].indice==right[iR].indice){
					update_connected(left[iL]);
				}
				//~ while(left[++iL].kmmer<=kR.kmmer ){}
				while(true){
					++iL;
					if(iL>=left.size()){break;}
					if(not (left[iL].kmmer<=kR.kmmer)){break;}
				}
			}
			if(right[iR].kmmer==kL.kmmer){
				go=false;
				if(not left[iL].indice==right[iR].indice){
					update_connected(right[iR]);
				}
				//~ while(right[++iR].kmmer<=kL.kmmer ){}
				while(true){
					++iR;
					if(iR>=right.size()){break;}
					if(not (right[iR].kmmer<=kL.kmmer)){break;}
				}
			}
			if(go){
				compaction(kL.indice,kR.indice,kL.kmmer);
				//~ to_compact.push_back(std::make_pair(kL.indice,kR.indice));
				}

		}else{
			if(kL.kmmer<kR.kmmer){
				//~ while(left[++iL].kmmer<kR.kmmer){}
				while(true){
					++iL;
					if(iL>=left.size()){break;}
					if(not (left[iL].kmmer<kR.kmmer)){break;}
				}
			}else{
				//~ while(right[++iR].kmmer<kL.kmmer){}
				while(true){
					++iR;
					if(iR>=right.size()){break;}
					if(not (right[iR].kmmer<kL.kmmer)){break;}
				}
			}
		}
	}

    //~ for (auto p: to_compact)
    //~ {
        //~ compaction(std::get<0>(p),std::get<1>(p));
    //~ }
}


template<size_t span>
bool graph3<span>::output(uint i){
    if (isNumber(unitigs[i][0]))
        return false;

    // Rayan: 
    // I remember thinking that indexed_left/indexed_right/connected_left/connected_right could not be properly set unless we did compactions at the end of the function above, i.e. using the 'to_compact' variable and doing 'for (auto p: to_compact) compaction(..)'.
    // but now, in the current version of the code, the 'to_compact' scheme is commented out
    // (due to https://github.com/GATB/gatb-core/commit/a33d9da4628ea18479798720db4207feb61e445b#diff-cd5d56912653ed079fe11a9d7d5ffbbb)
    // in order to restore the functionality of indexed_left/indexed_right/connected_left/connected_right, i would need to sit down and think again about it for a while. 
    // for now, just to be safe, i'm disabling pre_tip_cleaning, which wasn't used anyway. i believe this is the only location where indexed_left/indexed_right/connected_left/connected_right are used 
    // (in principle it could be possible to avoid doing the pass LinkTigs after bcalm/bglue, but it would need bglue to properly update the links. didn't code that.
#if 0
    if (pre_tip_cleaning)
    {
        if (indexed_left[i] && indexed_right[i])
        {
            if ((connected_left[i] && (!connected_right[i])) ||
                (connected_right[i] && (!connected_left[i])))
            {

                if (unitigs[i].size() < 3*(k+1)) // the spades tip length convention, to be tuned
                {
                    nb_pretips++;
                    //std::cout << "filtering tip " << unitigs[i] << " indexing l/r " << indexed_left[i] << " " << indexed_right[i] << " connected l/r " << connected_left[i] << " " << connected_right[i] << std::endl;
                    return false;
                }
            }
        }
    }
#endif 

    //std::cout << "returning seq " << unitigs[i] << " indexing l/r " << indexed_left[i] << " " << indexed_right[i] << " connected l/r " << connected_left[i] << " " << connected_right[i] << std::endl;
    return true;
}


template<size_t span>
bool graph3<span>::clear(){delete [] unitigs; delete [] unitigs_abundances;
    /* // nah, not needed. it wasn't the cause of the memory fragmentation, because even when the graph isn't constructed it still happens
    left.clear(); right.clear(); left.shrink_to_fit(); right.shrink_to_fit(); */
    return true;}


template<size_t span>
uint graph3<span>::size(){return indiceUnitigs;};


/* this function inserts sequences into the structure
 * while the code uses the term "unitigs", initially these sequences are just kmers (but later they will be unitigs)
 * sequences and their abundances are stored in a plain list
 * the index consists of two lists: left and right
 * both indices store tuples of the form (sequence index, canonical kmer)
 * the left index corresponds to kmers that are seen at the left of input sequence in forward strand, or on the right of unitigs in reverse strand
 * the right index is, well, the other ones. useful schema:
 *
 *    l                    r
 *   --->                --->
 *   ------------------------ input sequence (possibly a k-mer)
 *   <---                <---
 *     r                   l
 */
template<size_t span>
void graph3<span>::addtuple(tuple<string,uint,uint,uint>& tuple){
    // input tuple: <unitigs string, left minimizer, right minimizer, abundance>
	unitigs[indiceUnitigs]=get<0>(tuple);
	unitigs_abundances[indiceUnitigs].push_back(get<3>(tuple));

    bool debug = false;
    string debug_kmer = "GTTTTTTAGATTCTGAGTGGAACGATGAATG";

	if(minimizer==(get<1>(tuple))){
        indexed_left.push_back(true);
		typename graph3<span>::kmerType kmer1(beg2int128(unitigs[indiceUnitigs]));
		typename graph3<span>::kmerType kmer2(rcb(kmer1));
		if(kmer1<=kmer2){
            if (debug) if (unitigs[indiceUnitigs].compare(debug_kmer) == 0) std::cout << "for that seq " << unitigs[indiceUnitigs] << ", left kmer1 is " << kmer1 << " index " << indiceUnitigs << std::endl;
			left.push_back(kmerIndiceT<span>{indiceUnitigs,kmer1, SEQ_LEFT});
		}
		if(kmer2<=kmer1){
            if (debug) if (unitigs[indiceUnitigs].compare(debug_kmer) == 0) std::cout << "for that seq " << unitigs[indiceUnitigs] << ", left kmer2 is " << kmer1 << " index " << indiceUnitigs << std::endl;
			right.push_back(kmerIndiceT<span>{indiceUnitigs,kmer2, SEQ_LEFT});

		}
        // TODO probably to handle kmers that are their own reerse compelment, do:
        //  if (kmer2 < kmer1) instead of the "else"
        //  but i didnt test it yet, was chasing another bug, so let's implement that later
	}
    else
        indexed_left.push_back(false);
	if(minimizer==get<2>(tuple)){
        indexed_right.push_back(true);
		typename graph3<span>::kmerType kmer1(end2int128(unitigs[indiceUnitigs]));
		typename graph3<span>::kmerType kmer2(rcb(kmer1));

		if(kmer1<=kmer2){
            if (debug) if (unitigs[indiceUnitigs].compare(debug_kmer) == 0) std::cout << "for that seq " << unitigs[indiceUnitigs] << ", right kmer1 is " << kmer1 << " index " << indiceUnitigs << std::endl;
			right.push_back(kmerIndiceT<span>{indiceUnitigs,kmer1, SEQ_RIGHT});
		}
		if(kmer2<=kmer1){
            if (debug) if (unitigs[indiceUnitigs].compare(debug_kmer) == 0) std::cout << "for that seq " << unitigs[indiceUnitigs] << ", right kmer2 is " << kmer2 << std::endl;
			left.push_back(kmerIndiceT<span>{indiceUnitigs,kmer2, SEQ_RIGHT});

		}
	}
    else
        indexed_right.push_back(false);
	++indiceUnitigs;
}



}}}}
