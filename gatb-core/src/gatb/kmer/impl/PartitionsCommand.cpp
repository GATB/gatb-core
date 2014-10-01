/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <gatb/kmer/impl/PartitionsCommand.hpp>

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
PartitionsCommand<span>:: PartitionsCommand (
        Bag<Count>* solidKmers,
        Iterable<Type>&   partition,
        IHistogram*    histogram,
        ISynchronizer* synchro,
        u_int64_t&     totalKmerNbRef,
        size_t         abundance,
        IteratorListener* progress,
		PartiInfo<5> * pInfo,
		int parti,
		size_t      nbCores,
		size_t      kmerSize


    )
        : _abundance(abundance),
          _solidKmers(solidKmers, 10*1000, synchro),
          _partition(partition),
          _histogram(histogram),
          _progress(progress, synchro),
          _totalKmerNb(0),
          _totalKmerNbRef(totalKmerNbRef),
		  _pInfo(pInfo),
		  _parti_num(parti),
		  _nbCores(nbCores),
	      _kmerSize(kmerSize) {};

template<size_t span>
PartitionsCommand<span>::~PartitionsCommand()  {
        __sync_fetch_and_add (&_totalKmerNbRef, _totalKmerNb);
    };

template<size_t span>
void PartitionsCommand<span>::insert (const Count& kmer)
    {
        u_int32_t max_couv  = 2147483646;

        _totalKmerNb++;

        /** We should update the abundance histogram*/
        _histogram.inc (kmer.abundance);

        /** We check that the current abundance is in the correct range. */
        if (kmer.abundance >= this->_abundance && kmer.abundance <= max_couv)  {  this->_solidKmers.insert (kmer);  }
    };

/********************************************************************************/
/** in this scheme we count k-mers inside a partition by a hash table */
template<size_t span>
PartitionsByHashCommand<span>:: PartitionsByHashCommand (
        Bag<Count>*      solidKmers,
        Iterable<Type>&  partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        size_t          abundance,
        IteratorListener* progress,
        u_int64_t       hashMemory,
		PartiInfo<5> * pInfo,
		int parti,
		size_t      nbCores,
		size_t      kmerSize
														

														 
    )
        : PartitionsCommand<span> (solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress,pInfo,parti,nbCores,kmerSize), _hashMemory(hashMemory)  {};

template<size_t span>
void PartitionsByHashCommand<span>:: execute ()
    {
        size_t count=0;

        /** We need a map for storing part of solid kmers. */
        OAHash<Type> hash (_hashMemory);

        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL(it);

//        for (it->first(); !it->isDone(); it->next())
//        {
//            hash.increment (it->item());
//
//            /** Some display. */
//            if (++count == 100000)  {  this->_progress.inc (count);  count=0; }
//        }

		
		
		//with decompactage
		//should pass the kmer model here
		//superk
		u_int8_t		nbK, rem ;
		uint64_t compactedK;
		int ks = this->_kmerSize;
		Type un = 1;
		Type kmerMask = (un << (ks*2)) - un;
		size_t shift = 2*(ks-1);
		
        for (it->first(); !it->isDone(); it->next()) {
			Type superk = it->item();
			it->next();
			Type seedk = it->item();
			
			
			compactedK =  superk.getVal();
			nbK = (compactedK >> 56) & 255; // 8 bits poids fort = cpt //todo for large k values
			rem = nbK;
			

			Type temp = seedk;
			Type rev_temp = revcomp(temp,ks);
			Type newnt ;
			Type mink;
			
			for (int ii=0; ii< nbK; ii++,rem--) {
				mink = std::min (rev_temp, temp);
				
				//insert elem here
				hash.increment (mink);
								
				if(rem < 2) break;
				newnt =  ( superk >> ( 2*(rem-2)) ) & 3 ;
				
				temp = ((temp << 2 ) |  newnt   ) & kmerMask;
				newnt =  Type(comp_NT[newnt.getVal()]) ;
				rev_temp = ((rev_temp >> 2 ) |  (newnt << shift) ) & kmerMask;
			}
			
			
			
			
		}
		
		
		
		
        /** We loop over the solid kmers map. */
        Iterator < Abundance<Type> >* itKmerAbundance = hash.iterator();
        LOCAL (itKmerAbundance);

        for (itKmerAbundance->first(); !itKmerAbundance->isDone(); itKmerAbundance->next())
        {
            /** We may add this kmer to the solid kmers bag. */
           this->insert ((Count&) itKmerAbundance->item());
        }
    };

	

template<size_t span>
class SuperKReader
{
	typedef typename Kmer<span>::Type  Type;

public:

	void operator() (Type& elem)
	{
		//reading elems by pairs
		if(_first)
		{
			_superk = elem;
			
			_first = false;
		}
		else{
			_seedk = elem;

			
			uint64_t compactedK;
			
			compactedK =  _superk.getVal();
			u_int8_t nbK = (compactedK >> _shift_val    ) & 255; // 8 bits poids fort = cpt
			u_int8_t rem = nbK;
			
			Type temp = _seedk;
			Type rev_temp = revcomp(temp,_kmerSize);
			Type newnt ;
			Type mink, prev_mink;
			uint64_t idx;
			
			
			bool prev_which =  (temp < rev_temp );
			int kx_size = -1; //next loop start at ii=0, first kmer will put it at 0
			Type radix_kxmer_forward =  (temp & _mask_radix) >> ((_kmerSize - 4)*2);
			Type  first_revk, kinsert,radix_kxmer;
			
			if(!prev_which) first_revk = rev_temp;
			
			u_int8_t rid;

			//printf("init loop prev_which %i  kx_size %i \n",prev_which,kx_size);
			for (int ii=0; ii< nbK; ii++,rem--) {

				bool which =  (temp < rev_temp );
				mink = which ? temp : rev_temp;
				
			//	printf(" loop %i  which %i  kx_size %i \n",ii,which,kx_size);

				
				if(which != prev_which || kx_size >= _kx) // kxmer_size = 1
				{
					//output kxmer size kx_size,radix_kxmer
					//kx mer is composed of superKp[ii-1] superKp[ii-2] .. superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
					
					if(prev_which)
					{
						radix_kxmer = radix_kxmer_forward;
						kinsert = prev_mink;
					}
					else // si revcomp, le radix du kxmer est le debut du dernier kmer
					{
						//previous mink
						radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
						kinsert = first_revk;
					}
					
					//record kxmer
					rid = radix_kxmer.getVal();
					idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
					
					
			
					_radix_kmers[kx_size][rid][ idx] = kinsert << ((4-kx_size)*2);
					
					radix_kxmer_forward =  (mink & _mask_radix) >> _shift_radix;
					kx_size =0;
					
					if(!which) first_revk = rev_temp;

				}
				else
				{
					kx_size++;
				}
			
				prev_which = which ;
				prev_mink = mink;
				
				
				
				if(rem < 2) break; //no more kmers in this superkmer, the last one has just been eaten
				newnt =  ( _superk >> ( 2*(rem-2)) ) & 3 ;
				
				temp = ((temp << 2 ) |  newnt   ) & _kmerMask;
				newnt =  Type(comp_NT[newnt.getVal()]) ;
				rev_temp = ((rev_temp >> 2 ) |  (newnt << _shift) ) & _kmerMask;
			}
			
			
			
			//record last kxmer prev_mink et monk ?
			if(prev_which)
			{
				radix_kxmer = radix_kxmer_forward;
				kinsert = prev_mink;
			}
			else // si revcomp, le radix du kxmer est le debut du dernier kmer
			{
				//previous mink
				radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
				kinsert = first_revk;
			}
			
			//record kxmer
			rid = radix_kxmer.getVal();
			idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
			
					
			_radix_kmers[kx_size][rid][ idx] = kinsert << ((4-kx_size)*2);
			
			
			

			
			_first = true;

		}
		
	}
	
	SuperKReader (size_t kmerSize,  uint64_t * r_idx, std::vector <vector < vector<Type> > > & radix_kmers)
	: _first(true) ,_kmerSize (kmerSize), _r_idx (r_idx), _radix_kmers(radix_kmers), _kx(4)
	//model(model), pass(currentPass), nbPass(nbPasses), nbPartitions(partition->size()), nbWrittenKmers(0), _repart_table (repart_table),_mm(minim_size),
	 {
		 
		 Type un = 1;
		 _kmerMask = (un << (kmerSize*2)) - un;
		 _mask_radix  = Type((int64_t) 255);
		 _mask_radix =  _mask_radix << ((_kmerSize - 4)*2);
		 _shift = 2*(kmerSize-1);
		 _shift_val = un.getSize() -8;
		 _shift_radix = ((kmerSize - 4)*2); // radix is 4 nt long
	}

	
	private :
	size_t _kmerSize;
	size_t _shift ;
	size_t _shift_val ;
	size_t _shift_radix ;
	int _kx;
	std::vector < vector <vector<Type> > > & _radix_kmers;

	uint64_t * _r_idx ;
	bool _first;
	Type _superk, _seedk;
	
	
	Type _radix, _mask_radix ;
	
	
	Type _kmerMask;
};
	
	
//#define TIMP
/********************************************************************************/
/** in this scheme we count k-mers in a partition by sorting a vector*/
template<size_t span>
PartitionsByVectorCommand<span>:: PartitionsByVectorCommand (
        Bag<Count>*  solidKmers,
        Iterable<Type>&    partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        size_t          abundance,
        IteratorListener* progress,
		PartiInfo<5> * pInfo,
		int parti,
		size_t      nbCores,
		size_t      kmerSize
    )
        : PartitionsCommand<span> (solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress,pInfo,parti,nbCores,kmerSize)
          {};

template<size_t span>
void PartitionsByVectorCommand<span>:: execute ()
    {
		size_t _kx = 4 ;
		
		Dispatcher* parallel_dispatcher = 	new Dispatcher (this->_nbCores);
		//SerialDispatcher* parallel_dispatcher = 	new SerialDispatcher ();

		#ifdef TIMP
		 TimeInfo t;
		t.start ("lecture");
#endif
        /** We get the length of the current partition file. */
        size_t partitionLen = this->_partition.getNbItems();

        /** We check that we got something. */
        if (partitionLen == 0)  {  return;  /*throw Exception ("DSK: no solid kmers found");*/  }

		
		uint64_t _sum_nbxmer =0;
		uint64_t * r_idx  = (uint64_t *) calloc((_kx+1)*256,sizeof(uint64_t));
		memset(r_idx, 0 , sizeof(uint64_t)*(_kx+1)*256);
		
		radix_kmers.resize(_kx+1); // sapce for k0mer ... kxmer
		for (int xx=0; xx< (_kx+1); xx ++)
			radix_kmers[xx].resize(256);


		
		for (int xx=0; xx< (_kx+1); xx ++)
		for(int ii=0;ii< 256; ii++)
		{
		//	if( this->_pInfo->getNbKmer(this->_parti_num,ii,xx) !=0 )
		//	  printf("should resize  xmer %i rad %i  :  %llu \n",xx,ii, this->_pInfo->getNbKmer(this->_parti_num,ii,xx));
			radix_kmers[xx][ii].resize(this->_pInfo->getNbKmer(this->_parti_num,ii,xx));
			_sum_nbxmer +=  this->_pInfo->getNbKmer(this->_parti_num,ii,xx);
		}

		printf("--------- fillsolid parti num %i   nb kxmer / nbkmers      %lli / %lli     %f   with %zu nbcores -------\n",this->_parti_num,
			   _sum_nbxmer, this->_pInfo->getNbKmer(this->_parti_num),    (double) _sum_nbxmer /  this->_pInfo->getNbKmer(this->_parti_num),this->_nbCores );
		
        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL (it);
		
	
		//printf("decompacting the super kmers  \n");

		parallel_dispatcher->iterate (it, SuperKReader<span>  (this->_kmerSize, r_idx, radix_kmers), 10000); //must be even , reading by pairs
	
		//printf("done decompacting the super kmers  \n");

				
		//serial mode
//		for(int ii=0;ii< 256; ii++)
//		{
//				if(radix_kmers[ii].size() > 0)
//					std::sort (radix_kmers[ii].begin (), radix_kmers[ii].end ());
//		}
		
		#ifdef TIMP
		t.stop ("lecture");
		t.start ("tri");
#endif
	//	printf("sorting the kmers  \n");

		
		//parall mode
		vector<ICommand*> cmds;
		
		int nwork = 256 / this->_nbCores;
				
		// mettre dans le  SortCommand le master radix_kmers et range a traiter
		
		for (int xx=0; xx< (_kx+1); xx ++)
		{
			cmds.clear();
			//fill cmd work vector
			for(int tid=0;tid< this->_nbCores; tid++)
			{
				int deb = 0 + tid * nwork;
				int fin = (tid+1) * nwork -1; // thread will do inclusive range [begin -- end ]
				if(tid== this->_nbCores-1) fin = 255;
				
				ICommand* cmd = 0;
				
				cmd = new SortCommand<span> (radix_kmers[xx],deb,fin);
				cmds.push_back (cmd);
				
			}
			
			parallel_dispatcher->dispatchCommands (cmds, 0);
			
		}

//



        /** We sort the vector. */
	#ifdef TIMP
		t.stop ("tri");
		t.start ("output solid");
#endif
		
	//	printf("start scanning kx mers  \n");

		vector< KxmerPointer<span>  * > vec_pointer;
		
		int nbkxpointers = 453; //6 for k1 mer, 27 for k2mer, 112 for k3mer  453 for k4mer
		vec_pointer.resize(nbkxpointers);
		int best_p;

		std::priority_queue< kxp, std::vector<kxp>,kxpcomp > pq;

	
		u_int32_t abundance = 0;
		Type previous_kmer ;

		//puis loop sur tous les radix 0-255
		
		
	//	printf("start scanning \n");
			

		
		//init the pointers to the 6 arrays
		//	printf("init pointers  \n");
		int pidx =0;
		
		////-------------k0 pointers-----------/////////
		vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[0],0,0,0,255,this->_kmerSize ); // vec, prefix size, kxsize , radix min, radix max ,ksize
		
		////-------------k1 pointers-----------/////////
		//prefix0
		vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[1],0,1,0,255,this->_kmerSize );
		int lowr = 0;
		int maxr = 63;
		//prefix1
		for(unsigned int ii=0; ii<4; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[1],1,1,lowr,maxr,this->_kmerSize );
			lowr += 64;
			maxr += 64;
		}
		
		////-------------k2 pointers-----------/////////
		//prefix0
		vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[2],0,2,0,255,this->_kmerSize );
		//prefix1
		lowr = 0; maxr = 63;
		for(unsigned int ii=0; ii<4; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[2],1,2,lowr,maxr,this->_kmerSize );
			lowr += 64;
			maxr += 64;
		}
		
		//prefix2
		lowr = 0; maxr = 15;
		for(unsigned int ii=0; ii<16; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[2],2,2,lowr,maxr,this->_kmerSize );
			lowr += 16;
			maxr += 16;
		}
		
		////-------------k3 pointers-----------/////////

		//prefix0
		vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[3],0,3,0,255,this->_kmerSize );
		//prefix1
		lowr = 0; maxr = 63;
		for(unsigned int ii=0; ii<4; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[3],1,3,lowr,maxr,this->_kmerSize );
			lowr += 64;
			maxr += 64;
		}
		
		//prefix2
		lowr = 0; maxr = 15;
		for(unsigned int ii=0; ii<16; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[3],2,3,lowr,maxr,this->_kmerSize );
			lowr += 16;
			maxr += 16;
		}
		
		//prefix3
		lowr = 0; maxr = 3;
		for(unsigned int ii=0; ii<64; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[3],3,3,lowr,maxr,this->_kmerSize );
			lowr += 4;
			maxr += 4;
		}
		
		
		////-------------k4 pointers-----------/////////

		//prefix0
		vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[4],0,4,0,255,this->_kmerSize );
		//prefix1
		lowr = 0; maxr = 63;
		for(unsigned int ii=0; ii<4; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[4],1,4,lowr,maxr,this->_kmerSize );
			lowr += 64;
			maxr += 64;
		}
		
		//prefix2
		lowr = 0; maxr = 15;
		for(unsigned int ii=0; ii<16; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[4],2,4,lowr,maxr,this->_kmerSize );
			lowr += 16;
			maxr += 16;
		}
		
		//prefix3
		lowr = 0; maxr = 3;
		for(unsigned int ii=0; ii<64; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[4],3,4,lowr,maxr,this->_kmerSize );
			lowr += 4;
			maxr += 4;
		}
		
		//prefix4
		lowr = 0; maxr = 0;
		for(unsigned int ii=0; ii<256; ii++)
		{
			vec_pointer[pidx++] = new KxmerPointer<span> (radix_kmers[4],4,4,lowr,maxr,this->_kmerSize );
			lowr += 1;
			maxr += 1;
		}
		
		
		

		
		
		//printf("fill pq \n");

		//fill the  priority queue with the first elems
		for(int ii=0; ii<nbkxpointers; ii++)
		{
			if(vec_pointer[ii]->next())
				pq.push(kxp(ii,vec_pointer[ii]->value()));
		}

		
		if(pq.size() != 0) // everything empty, no kmer at all
		{
			//get first pointer
			best_p = pq.top().first ; pq.pop();
			

			previous_kmer = vec_pointer[best_p]->value();
			abundance = 1;
			
			
			//merge-scan all 'virtual' arrays and output counts
			while(1)
			{
				
				//go forward in this array or in new array of reaches end of this one
				if(! vec_pointer[best_p]->next())
				{

					
					//reaches end of one array
					if(pq.size() == 0) break; //everything done
					
					//otherwise get new best
					best_p = pq.top().first ; pq.pop();

				}
				
				if(vec_pointer[best_p]->value()!=previous_kmer )
				{

					//if diff, changes to new array, get new min pointer
					pq.push(kxp(best_p,vec_pointer[best_p]->value())); //push new val of this pointer in pq, will be counted later
					
					best_p = pq.top().first ; pq.pop();

					//if new best is diff, this is the end of this kmer
					if(vec_pointer[best_p]->value()!=previous_kmer )
					{

						this->insert (Count (previous_kmer, abundance) );
						abundance     = 1;
						previous_kmer = vec_pointer[best_p]->value();
					}
					else
					{

						abundance++;
					}
					
				}
				else
				{

					abundance++;
				}
			}
			
			
			//last elem
			this->insert (Count (previous_kmer, abundance) );
			
		}
		


		
		
		for(unsigned int ii=0; ii<nbkxpointers; ii++)
		{
			delete vec_pointer[ii];
		}
		


		
#ifdef TIMP
		t.stop ("output solid");

		cout << "lecture: " << t.getEntryByKey("lecture") << "  "
		<< "tri: " << t.getEntryByKey("tri") << "  "
		<< "output solid: " << t.getEntryByKey("output solid") << endl;
#endif

		/** We update the progress bar. */
        this->_progress.inc (this->_pInfo->getNbKmer(this->_parti_num) ); // this->_pInfo->getNbKmer(this->_parti_num)  kmers.size()

		//return;

    };

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class PartitionsCommand <KSIZE_1>;
template class PartitionsCommand <KSIZE_2>;
template class PartitionsCommand <KSIZE_3>;
template class PartitionsCommand <KSIZE_4>;

template class PartitionsByHashCommand <KSIZE_1>;
template class PartitionsByHashCommand <KSIZE_2>;
template class PartitionsByHashCommand <KSIZE_3>;
template class PartitionsByHashCommand <KSIZE_4>;

template class PartitionsByVectorCommand <KSIZE_1>;
template class PartitionsByVectorCommand <KSIZE_2>;
template class PartitionsByVectorCommand <KSIZE_3>;
template class PartitionsByVectorCommand <KSIZE_4>;


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
