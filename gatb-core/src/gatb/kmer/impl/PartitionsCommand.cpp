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
        IteratorListener* progress
    )
        : _abundance(abundance),
          _solidKmers(solidKmers, 10*1000, synchro),
          _partition(partition),
          _histogram(histogram),
          _progress(progress, synchro),
          _totalKmerNb(0),
          _totalKmerNbRef(totalKmerNbRef) {};

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
        u_int64_t       hashMemory
    )
        : PartitionsCommand<span> (solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress), _hashMemory(hashMemory)  {};

template<size_t span>
void PartitionsByHashCommand<span>:: execute ()
    {
        size_t count=0;

        /** We need a map for storing part of solid kmers. */
        OAHash<Type> hash (_hashMemory);

        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL(it);

        for (it->first(); !it->isDone(); it->next())
        {
            hash.increment (it->item());

            /** Some display. */
            if (++count == 100000)  {  this->_progress.inc (count);  count=0; }
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
        IteratorListener* progress
    )
        : PartitionsCommand<span> (solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress)
          {};

template<size_t span>
void PartitionsByVectorCommand<span>:: execute ()
    {
        /** We get the length of the current partition file. */
        size_t partitionLen = this->_partition.getNbItems();

        /** We check that we got something. */
        if (partitionLen == 0)  {  return;  /*throw Exception ("DSK: no solid kmers found");*/  }

        /** We resize our vector that will be filled with the partition file content.
         * NOTE: we add an extra item and we will set it to the maximum kmer value. */
        kmers.resize (1 + partitionLen);

        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL (it);
        size_t idx = 0;
        for (it->first(); !it->isDone(); it->next(), idx++) { kmers[idx] = it->item(); }

        /** We set the extra item to a max value, so we are sure it will sorted at the last location.
         * This trick allows to avoid extra treatment after the loop that computes the kmers abundance. */
        kmers[partitionLen] = ~0;

        /** We sort the vector. */
        std::sort (kmers.begin (), kmers.end ());

        u_int32_t abundance = 0;
        Type previous_kmer = kmers.front();

        /** We loop over the sorted solid kmers. */
        for (typename vector<Type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
        {
            if (*itKmers == previous_kmer)  {   abundance++;  }
            else
            {
                this->insert (Count (previous_kmer, abundance) );

                abundance     = 1;
                previous_kmer = *itKmers;
            }
        }

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
