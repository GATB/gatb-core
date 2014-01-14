/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

/** \file GraphHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HELPERS_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HELPERS_HPP_

/********************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {

/********************************************************************************/

class GraphHelper
{
public:

    /** */
    GraphHelper (const Graph& graph) : _graph(graph) {}

    /** */
    const Graph& getGraph() const { return _graph; }

    // simple paths traversal
    // invariant: the input kmer has no in-branching.
    // returns:
    //      1 if a good extension is found
    //      0 if a deadend was reached
    //     -1 if out-branching was detected
    //     -2 if no out-branching but next kmer has in-branching
    int simplePathAvance (const Node& node, Direction dir) const;
    int simplePathAvance (const Node& node, Direction dir, Edge& output) const;
    int simplePathAvance (const Node& node, Direction dir, kmer::Nucleotide& nt) const;

    /** */
    template<typename T> Graph::Iterator<T> simplePathIterator (const Node& node, Direction dir) const;

private:

    /** */
    const Graph& _graph;

    /** */
    Graph::Iterator<Node> getSimplePathNodeIterator (const Node& node, Direction dir) const;
};

/********************************************************************************/

template<> inline Graph::Iterator<Node> GraphHelper::simplePathIterator (const Node& node, Direction dir) const
{  return getSimplePathNodeIterator (node, dir);  }

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HELPERS_HPP_ */
