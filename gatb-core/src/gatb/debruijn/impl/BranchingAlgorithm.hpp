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
#include <gatb/debruijn/impl/GraphUnitigs.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

/** \brief Computation of the branching nodes of a Graph
 *
 * This class implements an algorithm that looks for branching nodes in the provided graph.
 *
 * A node N is branching <=> successors(N).size()!=1 || predecessors(N).size()!=1
 *
 * All found branching nodes are put into a storage object.
 *
 * Actually, this class is mainly used in the debruijn::impl::Graph class as a fourth step for
 * the de Bruijn graph creation.
 */
template <size_t span=KMER_DEFAULT_SPAN, typename Node=Node_t<>, typename Edge=Edge_t<Node_t<> >, typename Graph_t=Graph>
class BranchingAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::Type           Type;
    typedef typename kmer::impl::Kmer<span>::Count          Count;

    /** Constructor.
     * \param[in] graph : graph from which we look for branching nodes
     * \param[in] storage : storage where the found branching nodes will be put
     * \param[in] kind : kind of branching algorithm
     * \param[in] nb_cores : number of cores to be used; 0 means all available cores
     * \param[in] options : extra options
     */
    BranchingAlgorithm (
        const Graph_t& graph,
        tools::storage::impl::Storage& storage,
        tools::misc::BranchingKind  kind,
        size_t                      nb_cores = 0,
        tools::misc::IProperties*   options  = 0
    );

    /** Constructor.
     * \param[in] storage : retrieve the branching nodes from this storage.
     */
    BranchingAlgorithm (tools::storage::impl::Storage& storage);

    /** Destructor. */
    ~BranchingAlgorithm ();

    /** Get an option parser for branching parameters. Dynamic allocation, so must be released when no more used.
     * \return an instance of IOptionsParser.
     */
    static tools::misc::IOptionsParser* getOptionsParser ();

    /** \copydoc tools::misc::impl::Algorithm::execute */
    void execute ();

    /** Get the branching nodes as a collection of Count objects, ie couples of [kmer,abundance]
     * \return the collection of Count object.
     */
    tools::collections::Collection<Count>* getBranchingCollection() { return _branchingCollection; }

private:

    const Graph_t* _graph;

    tools::storage::impl::Storage& _storage;

    tools::misc::BranchingKind  _kind;

    tools::collections::Collection<Count>* _branchingCollection;
    void setBranchingCollection (tools::collections::Collection<Count>* branchingCollection)  {  SP_SETATTR(branchingCollection); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_ */

