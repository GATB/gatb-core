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

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <gatb/debruijn/impl/Terminator.hpp>

#include <cassert>

using namespace std;

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

#define DEBUG(a)   //a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::BranchingTerminatorTemplate (const GraphTemplate<Node,Edge,GraphDataVariant>& graph)
    : TerminatorTemplate<Node,Edge,GraphDataVariant> (graph)
{
    /** We loop over the branching nodes. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Iterator<BranchingNode_t<Node> > itBranching = this->_graph.GraphTemplate<Node,Edge,GraphDataVariant>::template iteratorBranching();
    for (itBranching.first(); !itBranching.isDone(); itBranching.next())
    {
        /** We add the current branching node into the map. */
        branching_kmers.insert (itBranching.item().kmer);
    }

    /** We finalize the map. We don't sort because branching are already sorted in the graph file. */
    branching_kmers.finalize (false);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::BranchingTerminatorTemplate (const BranchingTerminatorTemplate& terminator)
	: TerminatorTemplate<Node,Edge,GraphDataVariant>(terminator._graph), branching_kmers (terminator.branching_kmers)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::~BranchingTerminatorTemplate()
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::is_branching (Node& node) const
{
    return this->_graph.isBranching (node);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
int BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::getDelta (Edge& edge) const
{
         if (edge.direction == DIR_OUTCOMING && edge.from.strand == kmer::STRAND_FORWARD)  { return 0; }
    else if (edge.direction == DIR_OUTCOMING && edge.from.strand == kmer::STRAND_REVCOMP)  { return 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == kmer::STRAND_FORWARD)  { return 4; }
    else if (edge.direction == DIR_INCOMING  && edge.from.strand == kmer::STRAND_REVCOMP)  { return 0; }
    else { return -1; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
void BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::mark (Edge& edge)
{
    // BranchingTerminator ignores non-branching kmers
    if (!is_indexed (edge.from))  {   return;  }

    Value val=0;
    branching_kmers.get (edge.from.kmer, val);

    int delta = getDelta (edge);
    if (delta >= 0)
    {
        // set a 1 at the right NT & strand position
        val |= 1 << (edge.nt + delta);

        branching_kmers.set (edge.from.kmer,val); //was insert for Hash16
    }

    assert (is_marked(edge) == true);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::is_marked (Edge& edge)  const
{
    Value val = 0;
    int is_present = branching_kmers.get (edge.from.kmer, val);

    if (!is_present)  {   return false;  }

    int extension_nucleotide_marked = 0;

    int delta = getDelta (edge);
    if (delta >= 0)
    {
        // set a 1 at the right NT & strand position
        extension_nucleotide_marked = (val>>(edge.nt+delta))&1;
    }

    return  extension_nucleotide_marked == 1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
void BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::mark (Node& node)
{
    bool could_mark = false;

    // if it is a branching kmer, mark it directly (it may have no branching neighbor)
    if (is_indexed(node))
    {
        Value val = 0;
        branching_kmers.get (node.kmer, val);
        branching_kmers.set (node.kmer, val|(1<<8));
        could_mark = true;
    }

    /** We loop the neighbors edges of the current node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Edge> neighbors = this->_graph.template neighborsEdge(node.kmer);

    /** We loop the branching neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& e = neighbors[i];

        if (is_indexed(e.to)==false)  { continue; }

        /** We mark this edge (reversed first, in order to have the branching as the 'from' node) */
        Edge rev_e = this->_graph.reverse(e);
        mark (rev_e);

        could_mark = true;
    }

    if (could_mark)  {   assert(is_marked(node) == true);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::is_marked (Node& node) const
{
    // if it is a branching kmer, read marking directly (it may have no branching neighbor)
    if (is_indexed(node))  {   return is_marked_branching(node);  }

    /** We loop the neighbors edges of the current node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Edge> neighbors = this->_graph.neighborsEdge (node.kmer);

    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& e = neighbors[i];

        if  (is_indexed(e.to)==false)  { continue; }

        Edge rev_e = this->_graph.reverse(e);
        if (is_marked (rev_e) == true)  {  return true;  }
    }

    return false;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::is_marked_branching (Node& node) const
{
    Value val = 0;
    branching_kmers.get (node.kmer, val);
    return (val&(1<<8)) != 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
void BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::reset()
{
    branching_kmers.clear();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::is_indexed (Node& node) const
{
    return branching_kmers.contains (node.kmer);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
void BranchingTerminatorTemplate<Node,Edge,GraphDataVariant>::dump ()
{
}


/***********************/

template <typename Node, typename Edge, typename GraphDataVariant>
bool MPHFTerminatorTemplate<Node,Edge,GraphDataVariant>::is_marked (Node& node) const
{
    int status = this->_graph.queryNodeState(node);
    return status & 1 == 1;
}

template <typename Node, typename Edge, typename GraphDataVariant>
void MPHFTerminatorTemplate<Node,Edge,GraphDataVariant>::mark (Node& node) 
{
    int state = this->_graph.queryNodeState(node);
    state |= 1;
    this->_graph.setNodeState(node, state);
}

template <typename Node, typename Edge, typename GraphDataVariant>
void MPHFTerminatorTemplate<Node,Edge,GraphDataVariant>::reset() 
{
    this->_graph.resetNodeState();
}

// instantiation of GATB legacy compatibility 
template class MPHFTerminatorTemplate<Node, Edge, GraphDataVariant>; 
template class BranchingTerminatorTemplate<Node, Edge, GraphDataVariant>; 
template class TerminatorTemplate<Node, Edge, GraphDataVariant>; 


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
