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

/** \file SortingCountAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Counting kmers from a set of sequences
 */

#ifndef _SORTING_COUNT_ALGORITHM_HPP_
#define _SORTING_COUNT_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <string>

#include <gatb/kmer/impl/PartitionsCommand.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** \brief Class performing the kmer counting (also known as 'DSK')
 *
 * This class does the real job of counting the kmers from a reads database.
 *
 * This is a template class whose template argument is the kind of integer used for
 * kmers (integers on 64 bits, 128 bits, etc...)
 *
 * We define some template instantiations of this SortingCountAlgorithm; such an instantiation
 * does the real job of kmers counting. By defining several instantiations, we allow
 * to choose dynamically the correct class according to the user choice for kmer size
 * (remember that initial Minia version had to be re-compiled for different kmer size).
 */
template<size_t span=KMER_DEFAULT_SPAN>
class SortingCountAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::Type           Type;
    typedef typename kmer::impl::Kmer<span>::Count          Count;

    /** Constructor.*/
    SortingCountAlgorithm ();

    /** Constructor.*/
    SortingCountAlgorithm (
        tools::storage::impl::Storage* storage,
        gatb::core::bank::IBank* bank,
        size_t              kmerSize,
        size_t              abundance,
        u_int32_t           max_memory     = 0,
        u_int64_t           max_disk_space = 0,
        size_t              nbCores        = 0,
        size_t              partitionType  = 0,
        const std::string&  prefix         = "tmp.",
        const std::string&  histogramUri   = "histogram",
        gatb::core::tools::misc::IProperties* options = 0
    );

    /** Constructor.*/
    SortingCountAlgorithm (tools::storage::impl::Storage& storage);

    /** Destructor */
    virtual ~SortingCountAlgorithm ();

    /** operator=.*/
    SortingCountAlgorithm& operator= (const SortingCountAlgorithm& s);

    /** Process the kmers counting. It is mainly composed of a loop over the passes, and for each pass
     * 1) we build the partition files then 2) we fill the solid kmers file from the partitions.
     */
    void  execute ();

    /** Get the iterable over the computed solid kmers.
     * \return the solid kmers iterable. */
    tools::collections::Collection<Count>* getSolidKmers ()  { return _solidKmers; }

private:

    /** Compute several values, in particular the number of passes and partitions. */
    void configure (gatb::core::bank::IBank* bank);

    /** Fill partition files (for a given pass) from a sequence iterator.
     * \param[in] pass  : current pass whose value is used for choosing the partition file
     * \param[in] itSeq : sequences iterator whose sequence are cut into kmers to be split.
     */
    void fillPartitions (size_t pass, gatb::core::tools::dp::Iterator<gatb::core::bank::Sequence>* itSeq);

    /** Fill the solid kmers bag from the partition files (one partition after another one).
     * \param[in] solidKmers : bag to put the solid kmers into.
     */
    void fillSolidKmers (gatb::core::tools::collections::Bag<Count>*  solidKmers);

    /** */
    std::vector <size_t> getNbCoresList ();

    /** */
    tools::storage::impl::Storage* _storage;

    /** */
    gatb::core::bank::IBank* _bank;
    void setBank (gatb::core::bank::IBank* bank)  { SP_SETATTR(bank); }

    /** */
    tools::storage::impl::CollectionNode<Count>* _solidKmers;
    void setSolidKmers (tools::storage::impl::CollectionNode<Count>* solidKmers)
    {  _solidKmers = solidKmers;  }

    /** Shortcuts for the user input parameters. . */
    size_t      _kmerSize;
    size_t      _abundance;
    size_t      _partitionType;
    size_t      _nbCores;

    std::string _prefix;
    std::string _histogramUri;

    gatb::core::tools::dp::IteratorListener* _progress;
    void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }

    /** Values computed for algorithm parameterization. In particular, we have one value for the number
     * of passes and one value for the number of partitions.
     * Such values are computed both:
     *      - from system resources (file system resources, memory resources)
     *      - user preferences (max disk space, max memory)
     */
    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;
    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;
    u_int32_t _current_pass;

    gatb::core::tools::misc::IHistogram* _histogram;
    void setHistogram (gatb::core::tools::misc::IHistogram* histogram)  { SP_SETATTR(histogram); }

    /** Partitions management. */
    tools::storage::impl::Storage* _partitionsStorage;
    void setPartitionsStorage (tools::storage::impl::Storage* partitionsStorage)
    {
        SP_SETATTR(partitionsStorage);
    }

    tools::storage::impl::Partition<Type>* _partitions;
    void setPartitions (tools::storage::impl::Partition<Type>* partitions)  {  SP_SETATTR(partitions);  }

    u_int64_t _totalKmerNb;

};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _SORTING_COUNT_ALGORITHM_HPP_ */

