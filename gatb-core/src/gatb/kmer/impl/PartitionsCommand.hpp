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

#ifndef _PARTITIONSCOMMAND__HPP_
#define _PARTITIONSCOMMAND__HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/PartiInfo.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/bank/impl/Banks.hpp>

#include <gatb/tools/misc/impl/Pool.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

#include <queue>
#include <limits>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** \brief Class that counts the number of occurrences of a kmer in several banks.
 *
 * This class manages a vector of occurrences for a kmer, each index of the vector
 * holding the occurrences number for a bank.
 *
 * It is used during the process of kmers counting and eases the work for counting
 * kmers per bank.
 *
 * Note that it also implements the method \ref isSolid that tells whether the currently
 * counted kmer is solid or not, according to the solidity kind.
 */
class SolidityCounter
{
public:

    /** We define the integer type for getting abundance values. */
    typedef u_int32_t Int;

    /** We also define the maximum value for this integer type. */
    static Int MAX() { return std::numeric_limits<Int>::max(); }

    /** Constructor.
     * \param[in] solidityKind : tells how a kmer is considered as solid.
     * \param[in] threshold : range [min,max] for abundances
     * \param[in] nbBanks : number of banks parsed during kmer counting.
     */
    SolidityCounter (tools::misc::KmerSolidityKind solidityKind, const std::pair<Int,Int>& threshold, size_t nbBanks=1)
        : _solidityKind(solidityKind), _threshold(threshold), _abundancePerBank(nbBanks)
    {
        /** By convention, if max<min, we use max=MAX. */
        if (_threshold.second < _threshold.first)  {  _threshold.second = MAX(); }
    }

    /** Get the number of banks.
     * \return the number of banks. */
    size_t size() const  { return _abundancePerBank.size(); }

    /** Initialization of the counting for the current kmer. This method should be called
     * when a kmer is seen for the first time.
     * \param[in] idxBank : bank index where the new current kmer has been found. */
    void init (size_t idxBank=0)
    {
        for (size_t k=0; k<_abundancePerBank.size(); k++)  { _abundancePerBank[k]=0; }
        _abundancePerBank [idxBank]= 1;
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank=0)  {  _abundancePerBank [idxBank] ++;  }

    /** Get the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank
     * \return the abundance of the current kmer for the given bank. */
    Int operator[] (size_t idxBank) const  { return _abundancePerBank[idxBank]; }

    /** Tells whether the current kmer is solid or not. The computation is done
     * according to the solidity kind provided at constructor.
     * \return true if the current kmer is solid, false otherwise. */
    bool isSolid () const
    {
        /** Optimization. */
        if (size()==1)  { return (*this)[0] >= _threshold.first &&  (*this)[0] <= _threshold.second; }

        Int m = 0;

        /** By default we compute the sum. */
        for (size_t i=0; i<size(); i++)  { m += (*this)[i]; }

        switch (_solidityKind)
        {
            case tools::misc::KMER_SOLIDITY_MIN:
                m = ~0;  for (size_t i=0; i<size(); i++)  { if ((*this)[i] < m)  { m = (*this)[i]; } }
                break;

            case tools::misc::KMER_SOLIDITY_MAX:
                m = 0;  for (size_t i=0; i<size(); i++)  { if ((*this)[i] > m)  { m = (*this)[i]; } }
                break;

            case tools::misc::KMER_SOLIDITY_SUM:
            case tools::misc::KMER_SOLIDITY_DEFAULT:
            default:
                /** The sum is already computed as initialization value of 'm'. */
                break;
        }

        return (m >= _threshold.first) && (m <= _threshold.second);
    }

    /** Compute the sum of abundances of the current kmer for all the banks.
     * \return the sum of abundances. */
    Int computeSum () const
    {
        /** Optimization. */
        if (size()==1)  { return (*this)[0]; }

        Int sum=0; for (size_t k=0; k<_abundancePerBank.size(); k++)  { sum+=_abundancePerBank[k]; }  return sum;
    }

private:
    tools::misc::KmerSolidityKind _solidityKind;
    std::pair<Int,Int>            _threshold;
    std::vector<Int>              _abundancePerBank;
};

/********************************************************************************/
template<size_t span>
class PartitionsCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    /** Constructor. */
    PartitionsCommand (
        gatb::core::tools::collections::Bag<Count>*         solidKmers,
        gatb::core::tools::collections::Iterable<Type>&     partition,
        gatb::core::tools::misc::IHistogram*                histogram,
        gatb::core::system::ISynchronizer*                  synchro,
        u_int64_t&                                          totalKmerNbRef,
        std::pair<size_t,size_t>                            abundance,
        gatb::core::tools::dp::IteratorListener*            progress,
        tools::misc::impl::TimeInfo&                        timeInfo,
        PartiInfo<5>&                                       pInfo,
		int                                                 parti,
		size_t                                              nbCores,
		size_t                                              kmerSize,
		gatb::core::tools::misc::impl::MemAllocator&        pool,
	    size_t                                              cacheSize
    );

    /** Destructor. */
    ~PartitionsCommand();

protected:
    std::pair<size_t,size_t>                                _abundance;
    gatb::core::tools::collections::impl::BagCache<Count>   _solidKmers;
    gatb::core::tools::collections::Iterable<Type>&         _partition;
    gatb::core::tools::misc::impl::HistogramCache           _histogram;
    gatb::core::tools::dp::IteratorListener*                _progress;
    u_int64_t                                               _totalKmerNb;
    u_int64_t&                                              _totalKmerNbRef;
	PartiInfo<5>&                                           _pInfo;
	int                                                     _parti_num;
    size_t                                                  _nbCores;
	size_t                                                  _kmerSize;
	gatb::core::tools::misc::impl::MemAllocator&            _pool;
	
    void insert (const Count& kmer);

    void insert (const Type& kmer, const SolidityCounter& count);

    tools::misc::impl::TimeInfo& _globalTimeInfo;
    tools::misc::impl::TimeInfo  _timeInfo;
};

/********************************************************************************/
/** */
template<size_t span>
class PartitionsByHashCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication => R1: I'm afraid it's not possible. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    /** Constructor. */
    PartitionsByHashCommand (
        gatb::core::tools::collections::Bag<Count>*         solidKmers,
        gatb::core::tools::collections::Iterable<Type>&     partition,
        gatb::core::tools::misc::IHistogram*                histogram,
        gatb::core::system::ISynchronizer*                  synchro,
        u_int64_t&                                          totalKmerNbRef,
        std::pair<size_t,size_t>                            abundance,
        gatb::core::tools::dp::IteratorListener*            progress,
        tools::misc::impl::TimeInfo&                        timeInfo,
        PartiInfo<5>&                                       pInfo,
        int                                                 parti,
        size_t                                              nbCores,
        size_t                                              kmerSize,
        gatb::core::tools::misc::impl::MemAllocator&        pool,
        size_t                                              cacheSize,
        u_int64_t                                           hashMemory
    );

    void execute ();

private:
    u_int64_t _hashMemory;
};
		
/********************************************************************************/
/** */
template<size_t span>
class PartitionsByVectorCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    static const size_t KX = 4 ;

private:
    //used for the priority queue
    typedef std::pair<int, Type> kxp; //id pointer in vec_pointer , value
    struct kxpcomp { bool operator() (kxp l,kxp r) { return ((r.second) < (l.second)); } } ;

public:
    /** Constructor. */
    PartitionsByVectorCommand (
            gatb::core::tools::collections::Bag<Count>*         solidKmers,
            gatb::core::tools::collections::Iterable<Type>&     partition,
            gatb::core::tools::misc::IHistogram*                histogram,
            gatb::core::system::ISynchronizer*                  synchro,
            u_int64_t&                                          totalKmerNbRef,
            std::pair<size_t,size_t>                            abundance,
            gatb::core::tools::dp::IteratorListener*            progress,
            tools::misc::impl::TimeInfo&                        timeInfo,
            PartiInfo<5>&                                       pInfo,
            int                                                 parti,
            size_t                                              nbCores,
            size_t                                              kmerSize,
            gatb::core::tools::misc::impl::MemAllocator&        pool,
            size_t                                              cacheSize,
            tools::misc::KmerSolidityKind                       solidityKind,
            std::vector<size_t>&                                offsets
    );

    /** Destructor. */
    ~PartitionsByVectorCommand ();

    void execute ();

private:

    Type**     _radix_kmers;
    u_int8_t** _bankIdMatrix;
	uint64_t*  _radix_sizes;
	uint64_t*  _r_idx;

    tools::dp::IDispatcher* _dispatcher;

	void executeRead   ();
    void executeSort   ();
    void executeDump   ();

    std::vector<size_t> _nbItemsPerBankPerPart;

    tools::misc::KmerSolidityKind _solidityKind;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
