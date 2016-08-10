/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2015 INRIA
 *   Authors: R.Chikhi, G.Rizk, D.Lavenier, E.Drezen
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

// this class takes care of removing nodes from a graph
// not immediately, but only when flush() is called.

#ifndef _GATB_GRAPH_NODESDELETER_HPP_
#define _GATB_GRAPH_NODESDELETER_HPP_

/********************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/system/impl/System.hpp>
#include <vector>
#include <set>
#include <string>

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/


template <typename Node, typename Edge, typename Graph>
class NodesDeleter
{

    public:
        uint64_t nbNodes;
        std::vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
        std::set<Node> setNodesToDelete; 
        Graph &  _graph;
        int _nbCores;
        bool _verbose;
        bool useList, onlyListMethod;
        unsigned long explicitLimit;
        system::ISynchronizer* synchro;

    NodesDeleter(Graph&  graph, uint64_t nbNodes, int nbCores, bool verbose=true) : nbNodes(nbNodes), _graph(graph), _nbCores(nbCores), _verbose(verbose)
    {
        nodesToDelete.resize(nbNodes); // number of graph nodes // (!) this will alloc 1 bit per kmer.
        for (unsigned long i = 0; i < nbNodes; i++)
            nodesToDelete[i] = false;

        /* use explicit set of nodes as long as we don't have more than 10 M nodes ok? 
         * else resort to bit array
         * 10M nodes, assuming 128 bytes per nodes (generous), is 1 gig.
         * for k=21 it's actually closer to 24 bytes per node.
         * so, WARNING: enabling this feature uses more memory. 
         * todo: someday, introduce a --fast parameter that uses a bit more memory but does optimizations like this
         */
        //std::cout << "sizeof node: " << sizeof(Node) << std::endl;
        useList = true;
        onlyListMethod = false;
          
        // compute a fair amount of nodes that can be kept of memory
        // (before, was only 10M)
        explicitLimit = std::max((uint64_t)(nbNodes / 100), (uint64_t)10000000); 
        // for bacteria it's 10M
        // for human it's roughly 20M
        // for spruce it's roughly 300M
  
        // set insertions aren't thread safe, so let's use a synchronizer
        synchro = system::impl::System::thread().newSynchronizer();
    }

    ~NodesDeleter ()  { delete synchro; }

    bool get(uint64_t index)
    {
        return nodesToDelete[index];
    }
    
    bool get(Node &node)
    {
        if (useList)
        {
            return (setNodesToDelete.find(node) != setNodesToDelete.end());
        }
        //else
        {
            unsigned long index =_graph.nodeMPHFIndex(node);
            return get(index);
        }
    }

    void markToDelete(Node &node)
    {
        if (!onlyListMethod)
        {
            unsigned long index =_graph.nodeMPHFIndex(node);
            nodesToDelete[index] = true;
        }

        if (useList)
        {
            synchro->lock();
            setNodesToDelete.insert(node);
            synchro->unlock();
        }

        if ((!onlyListMethod) && (setNodesToDelete.size() > explicitLimit))
            useList = false;

    }

   // TODO speed: tell graph whenever all the neighbors of a node will be deleted too, that way, don't need to update their adjacency! 
    void flush()
    {
        if (useList)
        {
            if (_verbose)
                std::cout << "NodesDeleter mem usage prior to flush: " << (setNodesToDelete.size() * sizeof(Node)) / 1024 / 1024 << " MB" << std::endl;
            // sequential nodes deletion (no need for parallel here, as deleteNode() needs to be atomic anyway
            for (typename std::set<Node>::iterator it = setNodesToDelete.begin(); it != setNodesToDelete.end(); it++)
            {
                Node node = *it; // remove this line when (if ever?) deleteNode is const Node&
                _graph.deleteNode(node);
            }
        }
        else
        {
            _graph.deleteNodesByIndex(nodesToDelete, _nbCores, synchro);
        }
    }

};

/********************************************************************************/

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif 

