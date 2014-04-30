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
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/bank/impl/Banks.hpp>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::kmer::impl;


/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
class PartitionsCommand : public ICommand, public system::SmartPointer
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    PartitionsCommand (
        Bag<Count>* solidKmers,
        Iterable<Type>&   partition,
        IHistogram*    histogram,
        ISynchronizer* synchro,
        u_int64_t&     totalKmerNbRef,
        size_t         abundance,
        IteratorListener* progress
    );

    ~PartitionsCommand();

protected:
    size_t              _abundance;
    BagCache<Count>      _solidKmers;
    Iterable<Type>&      _partition;
    HistogramCache      _histogram;
    ProgressSynchro     _progress;
    u_int64_t           _totalKmerNb;
    u_int64_t&          _totalKmerNbRef;

    void insert (const Count& kmer);
};

/********************************************************************************/
/** */
template<size_t span>
class PartitionsByHashCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;


    PartitionsByHashCommand (
        Bag<Count>*      solidKmers,
        Iterable<Type>&  partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        size_t          abundance,
        IteratorListener* progress,
        u_int64_t       hashMemory
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

    PartitionsByVectorCommand (
        Bag<Count>*  solidKmers,
        Iterable<Type>&    partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        size_t          abundance,
        IteratorListener* progress
    );

    void execute ();

private:

    vector<Type> kmers;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/


#endif
