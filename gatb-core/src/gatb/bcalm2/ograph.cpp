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
void graph3<span>::compaction(uint iL,  uint iR){
	if(iR!=iL){
		bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
		if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
		if(b1){return compaction(stoi(unitigs[iL]),iR);}
		if(b2){return compaction(iL,stoi(unitigs[iR]));}

		typename graph3<span>::kmerType beg1;//(beg2int128(unitigs[iL])); // that kind of initialization isn't supported in LargeInt.
        beg1.setVal(beg2int128(unitigs[iL])); 
		typename graph3<span>::kmerType end2;//(end2int128(unitigs[iR]));
        end2.setVal(end2int128(unitigs[iR])); 

		if(beg1==end2){
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
            
            compact_abundances(iR,iL);
			return;
		}

		typename graph3<span>::kmerType endrc2;//(beg2int128rc(unitigs[iR]));
        endrc2.setVal(beg2int128rc(unitigs[iR])); 
		if(beg1==endrc2){
			reverseinplace2(unitigs[iR]);
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
            
            compact_abundances(iR,iL,true,false);
			return;
		}

		typename graph3<span>::kmerType beg2;//(rcb(endrc2));
        beg2.setVal(rcb(endrc2)); 
		typename graph3<span>::kmerType end1;//(end2int128(unitigs[iL]));
        end1.setVal(end2int128(unitigs[iL]));
		if(end1==beg2){
			unitigs[iL]+=(unitigs[iR].substr(k));
			unitigs[iR]=to_string(iL);

            compact_abundances(iL,iR);
			return;
		}

		typename graph3<span>::kmerType begrc2;//(rcb(end2));
        begrc2.setVal(rcb(end2));
/*            std::cout << "a : " << rcb(end2).toString(31) << std::endl;
            std::cout << "a=b: " << begrc2.toString(31) << std::endl;*/ // manifestation of a bug when LargeInt constructors are removed

		if(end1==begrc2){
			unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
			unitigs[iR]=to_string(iL);

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
void graph3<span>::debruijn(){
	sort(left.begin(),left.end(),comparator<span>());
	sort(right.begin(),right.end(),comparator<span>());
	uint iL(0),iR(0),sizeLeft(left.size()),sizeRight(right.size());
    typename graph3<span>::kmerType minusone;
    minusone.setVal(-1);
	left.push_back({0,minusone});
	right.push_back({0,minusone});

	kmerIndiceT<span> kL,kR;
	while(iL!=sizeLeft and iR!=sizeRight){
		kL=left[iL];
		kR=right[iR];
		if(kL.kmmer==kR.kmmer){
			bool go(true);
			++iL;++iR;
				if(left[iL].kmmer==kL.kmmer){
					go=false;
					while(left[++iL].kmmer==kL.kmmer){}
				}
				if(right[iR].kmmer==kL.kmmer){
					go=false;
					while(right[++iR].kmmer==kR.kmmer){}
				}
			if(go){compaction(kL.indice,kR.indice);}
		}else{
			if(kL.kmmer<kR.kmmer){
				while(left[++iL].kmmer==kL.kmmer){}
			}else{
				while(right[++iR].kmmer==kR.kmmer){}
			}
		}
	}
}


template<size_t span>
bool graph3<span>::output(uint i){return !isNumber(unitigs[i][0]);}


template<size_t span>
bool graph3<span>::clear(){delete [] unitigs; delete [] unitigs_abundances; return true;}


template<size_t span>
uint graph3<span>::size(){return indiceUnitigs;};


template<size_t span>
void graph3<span>::addtuple(tuple<string,uint,uint,uint>& tuple){
	unitigs[indiceUnitigs]=move(get<0>(tuple));
	unitigs_abundances[indiceUnitigs].push_back(get<3>(tuple));
	if(minimizer==(get<1>(tuple))){
		typename graph3<span>::kmerType kmer1(beg2int128(unitigs[indiceUnitigs]));
		typename graph3<span>::kmerType kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			left.push_back(kmerIndiceT<span>{indiceUnitigs,kmer1});
		}else{
			right.push_back(kmerIndiceT<span>{indiceUnitigs,kmer2});
		}
	}
	if(minimizer==get<2>(tuple)){
		typename graph3<span>::kmerType kmer1(end2int128(unitigs[indiceUnitigs]));
		typename graph3<span>::kmerType kmer2(rcb(kmer1));

		if(kmer1<kmer2){
			right.push_back(kmerIndiceT<span>{indiceUnitigs,kmer1});
		}else{
			left.push_back(kmerIndiceT<span>{indiceUnitigs,kmer2});
		}
	}
	++indiceUnitigs;
}


// void compareUnitigs(const string& fileFa,const string& fileDot){
// 	uint a(0),b(0),c(0),d(0);
// 	unordered_set<string> setFa,setDot;
// 	ifstream streamFa(fileFa),streamDot(fileDot);
// 	string seq;
// 	getline(streamFa,seq);
// 	while (!streamFa.eof()) {
// 		getline(streamFa,seq,'>');
// 		seq=seq.substr(0,seq.size()-1);
// 		setFa.insert(seq);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		getline(streamFa,seq);
// 		++c;
// 	}
// 	cout<<1<<endl;
// 	while (!streamDot.eof()){
// 		getline(streamDot,seq);
// 		transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
// 		seq=seq.substr(0,seq.size()-1);
// 		setDot.insert(seq);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		++d;
// 	}
// 	cout<<2<<endl;
// 	for(auto it(setFa.begin());it!=setFa.end();++it){
// 		if(setDot.count(*it)==0){
// 			++a;
// 		}
// 	}
// 	cout<<3<<endl;
// 	for(auto it(setDot.begin());it!=setDot.end();++it){
// 		if(setFa.count(*it)==0){
// 			++a;
// 		}
// 	}
// 	cout<<a<<" "<<b<<endl;
// 	cout<<c<<" "<<d<<endl;
// }
//
//
// void compareKmers(const string& fileFa,const string& fileDot){
// 	uint k(31);
// 	string kmer;
// 	uint a(0),b(0),c(0),d(0);
// 	unordered_set<string> setFa,setDot;
// 	ifstream streamFa(fileFa),streamDot(fileDot);
// 	string seq,inter,nimp;
//
//
//
// 	// cout<<1<<endl;
// 	while (!streamFa.eof()) {
// 		getline(streamFa,nimp);
// 		// cout<<"nimp"<<nimp<<endl;
// 		getline(streamFa,seq);
// 		// cout<<"seq"<<seq<<endl;
// 		point:
// 		char c=streamFa.peek();
// 		if(c=='>'){
// 			point2:
// 			// seq=seq.substr(0,seq.size());
// 			// for(uint j(0);(j)<seq.size();++j){
// 			// 	if(seq[j]!='A' and seq[j]!='C' and seq[j]!='T' and seq[j]!='G'){
// 			// 		cout<<seq<<endl;
// 			// 		cout<<"lol"<<endl;
// 			// 		exit(0);
// 			// 	}
// 			// }
// 			for (uint i = 0; i+k <=seq.size(); ++i) {
// 				kmer=seq.substr(i,k);
// 				// cout<<kmer<<endl;
// 				kmer=getRepresent(kmer);
// 				// if(setDot.count(kmer)==0){
// 				// 	++a;
// 				// }
// 				setFa.insert(kmer);
// 			}
// 		}else{
// 			if(!streamFa.eof()){
// 				// cout<<"inter"<<endl;
// 				// cout<<seq<<endl;
// 				getline(streamFa,inter);
// 				// cout<<inter<<endl;
// 				seq+=inter;
// 				goto point;
// 			}else{
// 				// cout<<"lol2"<<endl;
// 				goto point2;
// 			}
// 		}
// 	}
// 	cout<<2<<endl;
//
// 	while (!streamDot.eof()){
// 		getline(streamDot,seq);
// 		seq=seq.substr(0,k);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		if(setFa.count(getRepresent(seq))==0){
// 			cout<<seq<<endl;
// 			++a;
// 		}
// 	}
//
// 	// while (!streamDot.eof()){
// 	// 	getline(streamDot,seq);
// 	// 	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
// 	// 	seq=seq.substr(0,seq.size()-1);
// 	// 	// cout<<seq<<endl;
// 	// 	for (uint i = 0; i+k <=seq.size(); ++i) {
// 	// 		kmer=seq.substr(i,k);
// 	// 		// cout<<kmer<<endl;
// 	// 		kmer=getRepresent(kmer);
// 	// 		// setDot.insert(kmer);
// 	// 		if(setFa.count(kmer)==0){
// 	// 			++b;
// 	// 		}
// 	// 	}
// 	// 	// cout<<seq<<endl;
// 	// 	// cin.get();
// 	// 	// ++d;
// 	// }
// 	// for(auto it(setFa.begin());it!=setFa.end();++it){
// 	// 	if(setDot.count(*it)==0){
// 	// 		++a;
// 	// 	}
// 	// }
// 	cout<<3<<endl;
// 	// for(auto it(setDot.begin());it!=setDot.end();++it){
// 	// 	if(setFa.count(*it)==0){
// 	// 		++b;
// 	// 	}
// 	// }
// 	cout<<a<<" "<<b<<endl;
// 	cout<<c<<" "<<d<<endl;
// }

}}}}
