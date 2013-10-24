/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _SORTING_COUNT_ALGORITHM_HPP_
#define _SORTING_COUNT_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/Product.hpp>
#include <string>

#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** \brief Class performing the kmer counting.
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
template<typename ProductFactory, typename T> class SortingCountAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Constructor.*/
    SortingCountAlgorithm ();

    /** Constructor.*/
    SortingCountAlgorithm (
        tools::collections::impl::Product<ProductFactory>* product,
        gatb::core::bank::IBank* bank,
        size_t              kmerSize,
        size_t              nks,
        u_int32_t           max_memory     = 1000,
        u_int64_t           max_disk_space = 0,
        size_t              nbCores        = 0,
        size_t              partitionType  = 0,
        const std::string&  prefix         = "tmp.",
        const std::string&  histogramUri   = "histogram",
        gatb::core::tools::misc::IProperties* options = 0
    );

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
    tools::collections::Collection<Kmer<T> >* getSolidKmers ()  { return _solidKmers; }

private:

    typedef tools::collections::impl::ProductFileFactory  PartitionFactory;

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
    void fillSolidKmers (gatb::core::tools::collections::Bag<Kmer<T> >*  solidKmers);

    /** */
    std::vector <size_t> getNbCoresList ();

    /** */
    tools::collections::impl::Product<ProductFactory>* _product;

    /** */
    gatb::core::bank::IBank* _bank;
    void setBank (gatb::core::bank::IBank* bank)  { SP_SETATTR(bank); }

    /** */
    tools::collections::impl::CollectionNode<Kmer<T> >* _solidKmers;
    void setSolidKmers (tools::collections::impl::CollectionNode<Kmer<T> >* solidKmers)
    {  _solidKmers = solidKmers;  }

    /** Shortcuts for the user input parameters. . */
    size_t      _kmerSize;
    size_t      _nks;
    size_t      _partitionType;
    size_t      _nbCores;

    std::string _prefix;
    std::string _histogramUri;

    gatb::core::tools::dp::IteratorListener* _progress;
    void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }
    gatb::core::tools::dp::IteratorListener* getProgress ()  { return _progress; }

    /** */
    size_t getNks ()  { return _nks; }

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
    tools::collections::impl::Product<PartitionFactory>* _partitionsProduct;
    void setPartitionsProduct (tools::collections::impl::Product<PartitionFactory>* partitionsProduct)
    {
        SP_SETATTR(partitionsProduct);
    }

    tools::collections::impl::Partition<PartitionFactory, T>* _partitions;
    void setPartitions (tools::collections::impl::Partition<PartitionFactory, T>* partitions)  {  SP_SETATTR(partitions);  }

    u_int64_t _totalKmerNb;

    template<typename T1, typename T2> friend class PartitionsCommand;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _SORTING_COUNT_ALGORITHM_HPP_ */

