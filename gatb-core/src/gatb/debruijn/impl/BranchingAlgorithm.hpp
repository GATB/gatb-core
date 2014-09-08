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

/** \file BranchingAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Algorithm that computes the branching nodes of a De Bruijn graph
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

template <size_t span=KMER_DEFAULT_SPAN>
class BranchingAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::Type           Type;
    typedef typename kmer::impl::Kmer<span>::Count          Count;

    /** Constructor. */
    BranchingAlgorithm (
        const Graph& graph,
        tools::storage::impl::Storage& storage,
        tools::misc::BranchingKind  kind,
        size_t                      nb_cores = 0,
        tools::misc::IProperties*   options  = 0
    );

    /** Constructor. */
    BranchingAlgorithm (tools::storage::impl::Storage& storage);

    /** Destructor. */
    ~BranchingAlgorithm ();

    /** */
    void execute ();

    /** */
    tools::collections::Collection<Count>* getBranchingCollection() { return _branchingCollection; }

private:

    const Graph* _graph;

    tools::storage::impl::Storage& _storage;

    tools::misc::BranchingKind  _kind;

    tools::collections::Collection<Count>* _branchingCollection;
    void setBranchingCollection (tools::collections::Collection<Count>* branchingCollection)  {  SP_SETATTR(branchingCollection); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_ */

