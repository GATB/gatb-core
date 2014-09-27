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
		PartiInfo * pInfo,
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
		PartiInfo * pInfo,
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
		PartiInfo * pInfo,
		int parti,
		size_t      nbCores,
		size_t      kmerSize
    )
        : PartitionsCommand<span> (solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress,pInfo,parti,nbCores,kmerSize)
          {};

template<size_t span>
void PartitionsByVectorCommand<span>:: execute ()
    {
		
//		 TimeInfo t;
//		t.start ("lecture");
		
		//printf("fillsolid parti num %i \n",this->_parti_num);
        /** We get the length of the current partition file. */
        size_t partitionLen = this->_partition.getNbItems();

        /** We check that we got something. */
        if (partitionLen == 0)  {  return;  /*throw Exception ("DSK: no solid kmers found");*/  }

        /** We resize our vector that will be filled with the partition file content.
         * NOTE: we add an extra item and we will set it to the maximum kmer value. */
        //kmers.resize (1 + this->_pInfo->getNbKmer(this->_parti_num)); //ou reserve puis push back ?  le resize va appeler construc  de type pour rien ?
		// kmers.reserve (1 + this->_pInfo->getNbKmer(this->_parti_num));
		//kmers.reserve(4*partitionLen);
		

		radix_kmers.resize(256);
		for(int ii=0;ii< 256; ii++)
		{
			//printf("get %i %i \n", this->_parti_num,ii);
			radix_kmers[ii].reserve(1 + this->_pInfo->getNbKmer(this->_parti_num,ii));
		}
		
		Type mask_radix ((int64_t) 255);
		mask_radix =  mask_radix << ((this->_kmerSize - 4)*2); //get first 4 nt  of the kmers (heavy weight)
		Type radix ;

		
		//repartir les kmers ds les differents bins, faire tri de chaque ,
		//puis lecture de chaque separement possible aussi
		//puis version parall
		// puis version kx mer
		//puis gestion depass ram : multi pass  (mais comment estim ram max )? et cest fini
		
		
		//hmm do not know exactly how much will be needed, hence max mem can be  in worst case 2 times greater than necessary
		//todo save somewhere nb of elems per parti, to be able to resize accordingly  (then  kmers[idx] = ) instead of reserve
		
        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL (it);
        size_t idx = 0;
		
		//should pass the kmer model down here
		//superk
		u_int8_t		nbK, rem ;
		uint64_t compactedK;
		int ks = this->_kmerSize;

		Type un = 1;
		Type kmerMask = (un << (ks*2)) - un;
		
		//without decompactage
		
//		  for (it->first(); !it->isDone(); it->next(), idx++) {
//			kmers[idx] = it->item();
//		  }
		
		size_t shift = 2*(ks-1);
		
        for (it->first(); !it->isDone(); it->next()) {
		
			//here expand superk
			//kmers[idx] = it->item();
		
			
			// lecture 2 par 2
			Type superk = it->item();
			//printf("%s      %llx \n",superk.toString(ks).c_str(),superk.getVal());
			//
			it->next();
			
			Type seedk = it->item();
		//	printf("%s      %llx \n",seedk.toString(ks).c_str(),seedk.getVal());
			
			
			
			
			compactedK =  superk.getVal();
			nbK = (compactedK >> 56) & 255; // 8 bits poids fort = cpt //todo for large k values
			rem = nbK;
			
		//	printf("read new super k  %i  : \n",nbK);
			
			//temp =  *((uint64_t *)(&_seed) );
			Type temp = seedk;
			Type rev_temp = revcomp(temp,ks);
			Type newnt ;
			Type mink;
			
			for (int ii=0; ii< nbK; ii++,rem--) {
				
			//	revc = revcomp(temp,ks); //non, il faut calc plus intelligemment que ca, calc
				//revcomp au debut puis  par iteration comme temp plus bas

				mink = std::min (rev_temp, temp);
				
				radix =  (mink & mask_radix) >> ((ks - 4)*2);
				radix_kmers[radix.getVal()].push_back(mink);
				
				
				
				//kmers[idx] = mink; idx++;
				//kmers.push_back(mink);
				
			//  printf("%s   (%i / %i)   ( %lli )   (revc %lli   temp %lli )\n",mink.toString(ks).c_str(),rem,nbK,mink.getVal(),revc.getVal(),temp.getVal() );

				if(rem < 2) break;
				newnt =  ( superk >> ( 2*(rem-2)) ) & 3 ;
				
				temp = ((temp << 2 ) |  newnt   ) & kmerMask;
				newnt =  Type(comp_NT[newnt.getVal()]) ;
				rev_temp = ((rev_temp >> 2 ) |  (newnt << shift) ) & kmerMask;
			}
			
			
		}

		
        /** We set the extra item to a max value, so we are sure it will sorted at the last location.
         * This trick allows to avoid extra treatment after the loop that computes the kmers abundance. */
		//kmers[idx] = ~0; //GR : est on sur que ~0 va se convertir en k max qqsoit taille de kmer ?  et pas juste init avec max 64 bit int ?  (< max kmer)
		//et meme, ya pas un  prob si le max kmer arrive en vrai ?  kmer GGGGG..GGGG sauf si k 32 fait avec 128 bit
        //kmers[partitionLen] = ~0;
		///kmers.push_back(~0);
		
//		for(int ii=0;ii< 256; ii++)
//		{
//			if(radix_kmers[ii].size() > 0)
//				radix_kmers[ii].push_back(~0);
//		}
		
		
		
		//serial mode
//		for(int ii=0;ii< 256; ii++)
//		{
//				if(radix_kmers[ii].size() > 0)
//					std::sort (radix_kmers[ii].begin (), radix_kmers[ii].end ());
//		}
		
		//parall mode
		vector<ICommand*> cmds;
		
		int nwork = 256 / this->_nbCores;
				
		// mettre dans le  SortCommand le master radix_kmers et range a traiter
		for(int tid=0;tid< this->_nbCores; tid++)
		{
			int deb = 0 + tid * nwork;
			int fin = (tid+1) * nwork -1; // thread will do inclusive range [begin -- end ]
			if(tid== this->_nbCores-1) fin = 255;
			
			ICommand* cmd = 0;
			
			cmd = new SortCommand<span> (radix_kmers,deb,fin);
			cmds.push_back (cmd);
			
			
		}
		
		Dispatcher* disp = 	new Dispatcher ();
		disp->dispatchCommands (cmds, 0);

//
//		t.stop ("lecture");
//		t.start ("tri");

		
        /** We sort the vector. */
  //      std::sort (kmers.begin (), kmers.end ());
		
//		t.stop ("tri");
//		t.start ("output solid");

		//with radix bins
		for(int ii=0;ii< 256; ii++)
		{
			
			if(radix_kmers[ii].size() == 0) continue;
			
			u_int32_t abundance = 0;
			//Type previous_kmer = kmers.front();
			Type previous_kmer = radix_kmers[ii].front();
			
			/** We loop over the sorted solid kmers. */
			// for (typename vector<Type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
			for (typename vector<Type>::iterator itKmers = radix_kmers[ii].begin(); itKmers != radix_kmers[ii].end(); ++itKmers)
				
			{
				if (*itKmers == previous_kmer)  {   abundance++;  }
				else
				{
					this->insert (Count (previous_kmer, abundance) );
					
					abundance     = 1;
					previous_kmer = *itKmers;
				}
			}
			this->insert (Count (previous_kmer, abundance) ); //last elem // je sais pas pquoi avec radix il faut cette version le trick  insert max marche pas
			
			
		}

//		t.stop ("output solid");
//
//		cout << "lecture: " << t.getEntryByKey("lecture") << "  "
//		<< "tri: " << t.getEntryByKey("tri") << "  "
//		<< "output solid: " << t.getEntryByKey("output solid") << endl;
		
        /** We update the progress bar. */
        this->_progress.inc (kmers.size());
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
