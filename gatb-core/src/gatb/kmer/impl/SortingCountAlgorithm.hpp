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
#include <gatb/kmer/impl/BankKmers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/IterableHelpers.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <string>

#include <gatb/kmer/impl/PartitionsCommand.hpp>
#include <gatb/kmer/impl/LinearCounter.hpp>

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
 *
 * Actually, this class is mainly used in the debruijn::impl::Graph class as a first step for
 * the de Bruijn graph creation.
 */
template<size_t span=KMER_DEFAULT_SPAN>
class SortingCountAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
private:

	/** We define here a 'maximizer' in the mmers of a specific kmer. */
	struct  CustomMinimizer
	{
		template<class Model>  void init (const Model& model, typename Kmer<span>::Type& optimum) const
		{
			optimum = model.getKmerMax();
		}
		
		bool operator() (const typename Kmer<span>::Type& current, const typename Kmer<span>::Type& optimum) const
		{
			u_int64_t a = current.getVal() ;
			u_int64_t b = optimum.getVal() ;

			// test 3 consecutive identical nt
			//      u_int64_t a1 = a >>2 ;
			//      u_int64_t a2 = a >>4 ;
			//      a1 = (~( a ^ a1)) &  (~ (a ^ a2)) ;
			//      a1 =  ((a1 >>1) & a1) & _mask_0101 & _mmask_m1 ;
			//      if(a1 != 0) return false;

			// test 2 consecutive aa anywhere except beginning
			//      u_int64_t a1 =  ( ~a )  & (  ~a >>2 );
			
			// test si AA consecutif sauf au debut
			u_int64_t a1 =   ~(( a )   | (  a >>2 ));
			a1 = ((a1 >>1) & a1) & _mask_ma1 ;
			if (a1 != 0) return false;
			return (a<b);

			// return (current<optimum);
		}
		
		int        _mm;
		u_int64_t  _mmask_m1;
		u_int64_t  _mask_0101;
		u_int64_t  _mask_ma1;
		
		CustomMinimizer(int minim_size)
		{
			_mm        = minim_size;
			_mmask_m1  = (1 << ((_mm-2)*2)) -1 ;
			_mask_0101 = 0x5555555555555555  ;
			_mask_ma1  = _mask_0101 & _mmask_m1;
		}
		
		CustomMinimizer(const CustomMinimizer& cm)
		{
			_mm        = cm._mm;
			_mmask_m1  = cm._mmask_m1;
			_mask_0101 = cm._mask_0101;
			_mask_ma1  = cm._mask_ma1;
		}
	};
	
public:

	/* Shortcuts. */
	typedef typename Kmer<span>::ModelCanonical  ModelCanonical;
	typedef typename Kmer<span>::ModelDirect     ModelDirect;
	typedef typename kmer::impl::Kmer<span>::template ModelMinimizer <ModelCanonical> 	Model;
    typedef typename kmer::impl::Kmer<span>::Type  Type;
    typedef typename kmer::impl::Kmer<span>::Count Count;

    /** Constructor (default). */
    SortingCountAlgorithm ();

    /** Constructor.
     * \param[in] storage : storage where the solid kmers are saved
     * \param[in] bank : input bank from which solid kmers are counted
     * \param[in] kmerSize : size of kmers
     * \param[in] abundance : range [min,max] of abundances for the solidity criteria
     * \param[in] max_memory : max memory to be used (in MBytes)
     * \param[in] max_disk_space : max disk to be used (in MBytes)
     * \param[in] nbCores : number of cores to be used; 0 means all available cores
     * \param[in] solidityKind : criteria for solidity computation
     * \param[in] histogramMax : max number of values for the kmers histogram
     * \param[in] partitionType : kind of temporary partitions
     * \param[in] minimizerType : kind of minimizer (kmc2 or other)
     * \param[in] repartitionType : kind of repartition of minimizers into partitions
     * \param[in] minimizerSize : size of minimizers
     * \param[in] prefix : prefix used for output names
     * \param[in] options : extra options if any as a IProperties instance
     */
    SortingCountAlgorithm (
        tools::storage::impl::Storage* storage,
        gatb::core::bank::IBank* bank,
        size_t              kmerSize,
        std::pair<size_t,size_t> abundance,
        u_int32_t           max_memory     = 0,
        u_int64_t           max_disk_space = 0,
        size_t              nbCores        = 0,
        tools::misc::KmerSolidityKind solidityKind = tools::misc::KMER_SOLIDITY_DEFAULT,
        size_t              histogramMax   = 10000,
        size_t              partitionType  = 0,
        size_t              minimizerType  = 0,
        size_t              repartitionType= 0,
        size_t              minimizerSize  = 0,
        const std::string&  prefix         = "tmp.",
        gatb::core::tools::misc::IProperties* options = 0
    );

    /** Constructor.*/
    SortingCountAlgorithm (tools::storage::impl::Storage& storage);

    /** Destructor */
    virtual ~SortingCountAlgorithm ();

    /** operator=
     * \param[in] s : object to be copied. */
    SortingCountAlgorithm& operator= (const SortingCountAlgorithm& s);

    /** Get an option parser for kmers counting parameters. Dynamic allocation, so must be released when no more used.
     * \param[in] mandatory : tells whether an argument has to be mandatory
     * \return an instance of IOptionsParser. */
    static tools::misc::IOptionsParser* getOptionsParser (bool mandatory=true);

    /** Process the kmers counting. It is mainly composed of a loop over the passes, and for each pass
     * 1) we build the partition files then 2) we fill the solid kmers file from the partitions.
     */
    void  execute ();

    /** Get the iterable over the computed solid kmers.
     * \return the solid kmers iterable. */
    tools::storage::impl::Partition<Count>* getSolidCounts ()  { return _solidCounts; }

    /** Get the iterable over the computed solid kmers.
     * \return the solid kmers iterable. */
    tools::collections::Iterable<Type>* getSolidKmers   ()  { return _solidKmers;  }

    /** Return the storage group where dsk stuf is stored.
     * \return the Group instance.  */
    tools::storage::impl::Group& getStorageGroup() { return  (*_storage)("dsk"); }

    void setMinAutoThreshold(int value) { _min_auto_threshold = value; }

private:

    /** Compute several values, in particular the number of passes and partitions. */
    void configure (gatb::core::bank::IBank* bank);

    /** Fill partition files (for a given pass) from a sequence iterator.
     * \param[in] pass  : current pass whose value is used for choosing the partition file
     * \param[in] itSeq : sequences iterator whose sequence are cut into kmers to be split.
     */
    void fillPartitions (size_t pass, gatb::core::tools::dp::Iterator<gatb::core::bank::Sequence>* itSeq, PartiInfo<5>& pInfo);

    /** Fill the solid kmers bag from the partition files (one partition after another one).
     * \param[in] solidKmers : bag to put the solid kmers into.
     */
    void fillSolidKmers (PartiInfo<5>& pInfo);

    /** */
    std::vector <size_t> getNbCoresList (PartiInfo<5>& pInfo);

    /** */
    tools::storage::impl::Storage* _storage;

    /** */
    gatb::core::bank::IBank* _bank;
    void setBank (gatb::core::bank::IBank* bank)  { SP_SETATTR(bank); }

    /** */
    tools::collections::Iterable<Type>* _solidKmers;
    void setSolidKmers (tools::collections::Iterable<Type>* solidKmers)  { SP_SETATTR(solidKmers); }

    tools::storage::impl::Partition<Count>* _solidCounts;
    void setSolidCounts (tools::storage::impl::Partition<Count>* solidCounts)
    {
        SP_SETATTR(solidCounts);
        setSolidKmers (solidCounts ? new tools::collections::impl::IterableAdaptor<Count,Type,Count2TypeAdaptor> (*_solidCounts) : 0);
    }

    /** Shortcuts for the user input parameters. . */
    size_t      _kmerSize;
    std::pair<size_t,size_t> _abundance;
    size_t      _partitionType;
    size_t      _minimizerType;
    size_t      _repartitionType;
    size_t      _minimizerSize;
    size_t      _nbCores;
    size_t      _nbCores_per_partition;
    size_t      _nb_partitions_in_parallel;
    size_t      _minim_size;

    std::string _prefix;

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

    /** Get the memory size (in bytes) to be used by each item.
     * IMPORTANT : we may have to count both the size of Type and the size for the bank id. */
    int getSizeofPerItem () const { return Type::getSize()/8 + (_nbKmersPerPartitionPerBank.size()>1 ? sizeof(u_int8_t) : 0); }

    u_int64_t _totalKmerNb;

    struct Count2TypeAdaptor  {  Type& operator() (Count& c)  { return c.value; }  };

    tools::misc::impl::TimeInfo _fillTimeInfo;

    u_int64_t _estimatedDistinctKmerNb;
    bool _flagEstimateNbDistinctKmers; // whether we estimate the number of distinct kmers beforehand
	
    std::pair<int,int> _partCmdTypes;

    BankStats _bankStats;

    std::vector <std::vector<size_t> > _nbKmersPerPartitionPerBank;

    tools::misc::KmerSolidityKind _solidityKind;

    int _min_auto_threshold; // used for histogram.compute_threshold() : prevents the auto_cutoff from being below this value. Default =3
};
	
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _SORTING_COUNT_ALGORITHM_HPP_ */

