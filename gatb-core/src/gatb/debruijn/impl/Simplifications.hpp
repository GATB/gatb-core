/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
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

#ifndef _GATB_GRAPH_SIMPLIFICATION_HPP_
#define _GATB_GRAPH_SIMPLIFICATION_HPP_

/********************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>
#include <vector>
#include <set>
#include <string>


/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

/** */
template<typename Node, typename Edge, typename GraphDataVariant>
class Simplifications : public system::SmartPointer
{
public:

    Simplifications (/*const, removed because of cacheNonSimpleNodes calling a setStats */ GraphTemplate<Node,Edge,GraphDataVariant> & graph, int nbCores, bool verbose = false);

    void simplify(); // perform many rounds of all simplifications, as in Minia

    unsigned long removeTips();
    unsigned long removeBulges();
    unsigned long removeErroneousConnections();

    double getSimplePathCoverage(Node node, Direction dir, unsigned int* pathLen = NULL, unsigned int maxLength = 0);
    double getMeanAbundanceOfNeighbors(Node branchingNode, Node nodeToExclude);
    bool satisfyRCTC(std::vector<Node>& nodes, double RCTCcutoff);

    int _nbTipRemovalPasses;
    int _nbBubbleRemovalPasses;
    int _nbBulgeRemovalPasses;
    int _nbECRemovalPasses;
    
    std::string tipRemoval, bubbleRemoval, ECRemoval;
    

protected:
    /*const*/ GraphTemplate<Node,Edge,GraphDataVariant> &  _graph;
    int _nbCores;
    uint64_t nbNodes;
    uint64_t cutoffEvents;

    bool _firstNodeIteration;
    bool _verbose;

    std::string path2string(Direction dir, Path_t<Node> p, Node endNode);
    double path2abundance(Direction dir, Path_t<Node> p, Node endNode);

    void heuristic_most_covered_path(Direction dir, Node& startingNode, Node& endingNode, 
                                    int traversal_depth, int& success, double& mean_abundance,
                                    Path_t<Node> &res_path,
                                    unsigned int backtrackingLimit = 0, Node *avoidFirstNode = nullptr,
                                    bool most_covered = true, bool cached = false);
    void heuristic_most_covered_path(Direction dir, Node& startingNode, Node& endingNode, 
                                    int traversal_depth, Path_t<Node>& current_path, std::set<typename Node::Value>& usedNode, int& success, std::vector<int>& abundances,
                                    unsigned int backtrackingLimit, Node *avoidFirstNode, 
                                    bool most_covered, Path_t<Node> &res_path,
                                    unsigned long &nbCalls);
    void heuristic_most_covered_path_cached(Direction dir, Node& startingNode, Node& endNode, 
                                    int traversal_depth, Path_t<Node>& current_path, std::set<typename Node::Value>& usedNode, int& success, double& mean_abundance,
                                    unsigned int backtrackingLimit, Node *avoidFirstNode, 
                                    bool most_covered, Path_t<Node> &res_path,
                                    unsigned long &nbCalls);


    std::vector<bool> interestingNodes;
};

/********************************************************************************/

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif 

