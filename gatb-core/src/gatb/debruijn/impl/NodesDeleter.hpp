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


template <typename Node, typename Edge, typename GraphDataVariant>
class NodesDeleter
{

    public:
        uint64_t nbNodes;
        std::vector<bool> nodesToDelete; // don't delete while parallel traversal, do it afterwards
        std::set<Node> setNodesToDelete; 
        const GraphTemplate<Node,Edge,GraphDataVariant> &  _graph;
        int _nbCores;
        bool useList;
        unsigned long explicitLimit;
        system::ISynchronizer* synchro;

    NodesDeleter(const GraphTemplate<Node,Edge,GraphDataVariant> &  graph, uint64_t nbNodes, int nbCores) : nbNodes(nbNodes), _graph(graph), _nbCores(nbCores)
    {
        nodesToDelete.resize(nbNodes); // number of graph nodes // (!) this will alloc 1 bit per kmer.
        for (unsigned long i = 0; i < nbNodes; i++)
            nodesToDelete[i] = false;

        /* use explicit set of nodes as long as we don't have more than 10 M nodes ok? 
         * else resort to bit array
         * 10M nodes, assuming 128 bytes per nodes (generous), is 1 gig.
         * for k=21 it's 24 bytes.
         * so, WARNING: enabling this feature uses more memory. 
         * todo: someday, introduce a --fast parameter that uses a bit more memory but does optimizations like this
         */
        //std::cout << "sizeof node: " << sizeof(Node) << std::endl;
        useList = true;
        explicitLimit = 10000000; 
            
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
        unsigned long index =_graph.nodeMPHFIndex(node);
        nodesToDelete[index] = true;

        if (useList)
        {
            synchro->lock();
            setNodesToDelete.insert(node);
            synchro->unlock();
        }

        if (setNodesToDelete.size() > explicitLimit)
            useList = false;

    }
    
    void flush()
    {
        if (useList)
        {
            // sequential nodes deletion (no need for parallel here, as deleteNode() needs to be atomic anyway
            for (typename std::set<Node>::iterator it = setNodesToDelete.begin(); it != setNodesToDelete.end(); it++)
                _graph.deleteNode(*it);
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

