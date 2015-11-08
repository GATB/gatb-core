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

#include <gatb/debruijn/impl/Frontline.hpp>

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
// a frontline is a set of nodes having equal depth in the BFS
template <typename Node, typename Edge, typename GraphDataVariant>
FrontlineTemplate<Node,Edge,GraphDataVariant>::FrontlineTemplate (
    Direction         direction,
    const GraphTemplate<Node,Edge,GraphDataVariant>&      graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&       terminator,
    Node&       startingNode
) :
    _direction(direction), _graph(graph), _terminator(terminator), _depth(0),
    _all_involved_extensions(0)
{
    _already_frontlined.insert (startingNode.kmer);

    _frontline.push (NodeNt<Node>(startingNode, kmer::NUCL_UNKNOWN));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// a frontline is a set of nodes having equal depth in the BFS
template <typename Node, typename Edge, typename GraphDataVariant>
FrontlineTemplate<Node,Edge,GraphDataVariant>::FrontlineTemplate (
    Direction         direction,
    const GraphTemplate<Node,Edge,GraphDataVariant>&      graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&       terminator,
    Node&       startingNode,
    Node&       previousNode,
    std::set<Node>*   all_involved_extensions
) :
    _direction(direction), _graph(graph), _terminator(terminator), _depth(0),
    _all_involved_extensions(all_involved_extensions)
{
    _already_frontlined.insert (startingNode.kmer);
    _already_frontlined.insert (previousNode.kmer);

    _frontline.push (NodeNt<Node>(startingNode, kmer::NUCL_UNKNOWN));
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
bool FrontlineTemplate<Node,Edge,GraphDataVariant>::go_next_depth()
{
    // extend all nodes in this frontline simultaneously, creating a new frontline
    stopped_reason=NONE;
    queue_nodes new_frontline;

    while (!this->_frontline.empty())
    {
        /** We get the first item of the queue and remove it from the queue. */
        NodeNt<Node> current_node = this->_frontline.front();
        _frontline.pop();

        /** We check whether we use this node or not. we always use the first node at depth 0 */
        if (_depth > 0 && check(current_node.node) == false)  { return false; }

        /** We loop the neighbors edges of the current node. */
        typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Edge> edges = _graph.neighborsEdge (current_node.node, _direction);

        for (size_t i=0; i<edges.size(); i++)
        {
            /** Shortcuts. */
            Edge& edge     = edges[i];
            Node& neighbor = edge.to;

            // test if that node hasn't already been explored
            if (_already_frontlined.find (neighbor.kmer) != _already_frontlined.end())  { continue; }

            // if this bubble contains a marked (branching) kmer, stop everyone at once (to avoid redundancy)
            //if (_terminator.isEnabled() && _terminator.is_branching (neighbor) &&  _terminator.is_marked_branching(neighbor))   // legacy, before MPHFTerminator
            if (_terminator.isEnabled() && _terminator.is_marked(neighbor))   // to accomodate MPHFTerminator
            {  
                stopped_reason=FrontlineTemplate<Node,Edge,GraphDataVariant>::MARKED;
                return false;  
            }

            // propagate information where this node comes from
            kmer::Nucleotide from_nt = (current_node.nt == kmer::NUCL_UNKNOWN) ? edge.nt : current_node.nt;

            /** We add the new node to the new front line. */
            new_frontline.push (NodeNt<Node> (neighbor, from_nt));

            /** We memorize the new node. */
            _already_frontlined.insert (neighbor.kmer);

            // since this extension is validated, insert into the list of involved ones
            if (_all_involved_extensions != 0)  {  _all_involved_extensions->insert (neighbor);  }
        }
    }

    _frontline = new_frontline;
    ++_depth;

    return true;
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
FrontlineBranchingTemplate<Node,Edge,GraphDataVariant>::FrontlineBranchingTemplate (
    Direction         direction,
    const GraphTemplate<Node,Edge,GraphDataVariant>&      graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&       terminator,
    Node&       startingNode,
    Node&       previousNode,
    std::set<Node>*   all_involved_extensions
)  : FrontlineTemplate<Node,Edge,GraphDataVariant>(direction,graph,terminator,startingNode,previousNode,all_involved_extensions)
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
FrontlineBranchingTemplate<Node,Edge,GraphDataVariant>::FrontlineBranchingTemplate (
    Direction         direction,
    const GraphTemplate<Node,Edge,GraphDataVariant>&      graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&       terminator,
    Node&       startingNode
) : FrontlineTemplate<Node,Edge,GraphDataVariant>(direction,graph,terminator,startingNode)
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
// new code, not in monument, to detect any in-branching longer than 3k
template <typename Node, typename Edge, typename GraphDataVariant>
bool FrontlineBranchingTemplate<Node,Edge,GraphDataVariant>::check (Node& node)
{
	/** We reverse the node for the inbranching path. */
    Node actual = this->_graph.reverse(node);

    /** We loop the neighbors nodes of the current node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Node> neighbors = this->_graph.neighbors(actual, (this->_direction));

    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Node& neighbor = neighbors[i];

        // only check in-branching from kmers not already frontlined
        // which, for the first extension, includes the previously traversed kmer (previous_kmer)
        // btw due to avance() invariant, previous_kmer is always within a simple path
        if (this->_already_frontlined.find (neighbor.kmer) != this->_already_frontlined.end())  {   continue;  }

        // create a new frontline inside this frontline to check for large in-branching (i know, we need to go deeper, etc..)
        FrontlineTemplate<Node,Edge,GraphDataVariant> frontline (this->_direction, this->_graph, this->_terminator, neighbor, actual, this->_all_involved_extensions);

        do  {
            bool should_continue = frontline.go_next_depth();

            if (!should_continue)  
            {  
                this->stopped_reason = FrontlineTemplate<Node,Edge,GraphDataVariant>::IN_BRANCHING_OTHER;
                break;
            }

            // don't allow a depth > 3k
            if (frontline.depth() > 3 * this->_graph.getKmerSize())  
            {  
                this->stopped_reason = FrontlineTemplate<Node,Edge,GraphDataVariant>::IN_BRANCHING_DEPTH;
                break;
            }

            // don't allow a breadth too large
            if (frontline.size() > 10)  
            {  
                this->stopped_reason = FrontlineTemplate<Node,Edge,GraphDataVariant>::IN_BRANCHING_BREADTH;
                break;
            }

            // stopping condition: no more in-branching
            if (frontline.size() == 0)  {  break;  }
        }
        while (1);

        // found large in-branching
        if (frontline.size() > 0)  {  return false;  }
    }

    // didn't find any in-branching
    return true;
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
FrontlineReachableTemplate<Node,Edge,GraphDataVariant>::FrontlineReachableTemplate(
    Direction         direction,
    const GraphTemplate<Node,Edge,GraphDataVariant>&      graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&       terminator,
    Node&       startingNode,
    Node&       previousNode,
    std::set<Node>*   all_involved_extensions
)  : FrontlineTemplate<Node,Edge,GraphDataVariant> (direction,graph,terminator,startingNode,previousNode,all_involved_extensions)
{
}



template <typename Node, typename Edge, typename GraphDataVariant>
bool FrontlineReachableTemplate<Node,Edge,GraphDataVariant>::check (Node& node)
{
	/** We reverse the node for the inbranching path. */
    Node actual = this->_graph.reverse(node);

    /** neighbors nodes of the current node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Node> neighbors = this->_graph.neighbors(actual, (this->_direction));

    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Node& neighbor = neighbors[i];
        if (this->_already_frontlined.find (neighbor.kmer) == this->_already_frontlined.end())  {
            checkLater.insert(neighbor);
           //return false;   // strict
        }
    }
    return true;
}

template <typename Node, typename Edge, typename GraphDataVariant>
bool FrontlineReachableTemplate<Node,Edge,GraphDataVariant>::isReachable()
{
   for (typename std::set<Node>::iterator itNode = checkLater.begin(); itNode != checkLater.end(); itNode++)
   {
        if (this->_already_frontlined.find((*itNode).kmer) == this->_already_frontlined.end())
            return false;

   }
   return true;
}

// legacy GATB compatbility
template class FrontlineTemplate<Node, Edge, GraphDataVariant>; 
template class FrontlineBranchingTemplate<Node, Edge, GraphDataVariant>; 
template class FrontlineReachableTemplate<Node, Edge, GraphDataVariant>; 


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
