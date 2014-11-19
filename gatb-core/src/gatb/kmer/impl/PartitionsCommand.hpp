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

#include <queue>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

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
        size_t                                              abundance,
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
    size_t                                                  _abundance;
    gatb::core::tools::collections::impl::BagCache<Count>   _solidKmers;
    gatb::core::tools::collections::Iterable<Type>&         _partition;
    gatb::core::tools::misc::impl::HistogramCache           _histogram;
    gatb::core::tools::misc::impl::ProgressSynchro          _progress;
    u_int64_t                                               _totalKmerNb;
    u_int64_t&                                              _totalKmerNbRef;
	PartiInfo<5>&                                           _pInfo;
	int                                                     _parti_num;
    size_t                                                  _nbCores;
	size_t                                                  _kmerSize;
	gatb::core::tools::misc::impl::MemAllocator&            _pool;
	
    void insert (const Count& kmer);

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
        size_t                                              abundance,
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
            size_t                                              abundance,
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
    ~PartitionsByVectorCommand ();

    void execute ();

private:
	Type**    _radix_kmers;
	uint64_t* _radix_sizes;
	uint64_t* _r_idx;

    tools::dp::impl::Dispatcher* _dispatcher;

	void executeRead   ();
    void executeSort   ();
    void executeDump   ();
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
