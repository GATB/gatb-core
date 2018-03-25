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
template<typename GraphType, typename Node, typename Edge>
class Simplifications : public system::SmartPointer
{
public:

    Simplifications (/*const, removed because of cacheNonSimpleNodes calling a setStats */ GraphType * /* set as a pointer because could be null*/ graph, int nbCores, bool verbose = false);

    void simplify(); // perform many rounds of all simplifications, as in Minia

    unsigned long removeTips();
    unsigned long removeBulges();
    unsigned long removeErroneousConnections();

    double getSimplePathCoverage(Node node, Direction dir, unsigned int* pathLen = NULL, unsigned int maxLength = 0);
    double getMeanAbundanceOfNeighbors(Node& branchingNode, Node nodeToExclude);
    bool satisfyRCTC(double pathAbundance, Node& lastTipNode, double RCTCcutoff, Direction dir);

    int _nbTipRemovalPasses;
    int _nbBubbleRemovalPasses;
    int _nbBulgeRemovalPasses;
    int _nbECRemovalPasses;
    
    std::string tipRemoval, bubbleRemoval, ECRemoval;
    bool _doTipRemoval, _doBulgeRemoval, _doECRemoval;
   
    /* now exposing some parameters */
    double _tipLen_Topo_kMult;
    double _tipLen_RCTC_kMult;
    double _tipRCTCcutoff;

    double       _bulgeLen_kMult;
    unsigned int _bulgeLen_kAdd;
    unsigned int _bulgeAltPath_kAdd;
    unsigned int _bulgeAltPath_covMult;

    double _ecLen_kMult;
    double _ecRCTCcutoff;

protected:
    /*const*/ GraphType &  _graph;
    int _nbCores;
    uint64_t nbNodes;
    uint64_t cutoffEvents;

    bool _firstNodeIteration;
    bool _verbose;

    std::string path2string(Direction dir, Path_t<Node> p, Node endNode);
    double path2abundance(Direction dir, Path_t<Node> p, Node endNode, unsigned int skip_first = 0, unsigned int skip_last = 0);

    void heuristic_most_covered_path(Direction dir, Node& startingNode, Node& endingNode, 
                                    int traversal_depth, int& success, double& mean_abundance,
                                    Path_t<Node> &res_path,
                                    unsigned int backtrackingLimit = 0, Node *avoidFirstNode = NULL /*nullptr ideally, but want old gcc compatibility at least for headers*/,
                                    bool most_covered = true, bool kmer_version = false);
    // kmer version
    void heuristic_most_covered_path_old(Direction dir, Node& startingNode, Node& endingNode, 
                                    int traversal_depth, Path_t<Node>& current_path, std::set<Node>& usedNode, int& success, std::vector<int>& abundances,
                                    unsigned int backtrackingLimit, Node *avoidFirstNode, 
                                    bool most_covered, Path_t<Node> &res_path,
                                    unsigned long &nbCalls);
    // in-between version, towards unitigs
    void heuristic_most_covered_path(Direction dir, Node& startingNode, Node& endNode, 
                                    int traversal_depth, Path_t<Node>& current_path, std::set<Node>& usedNode, int& success, double& mean_abundance,
                                    unsigned int backtrackingLimit, Node *avoidFirstNode, 
                                    bool most_covered, Path_t<Node> &res_path,
                                    unsigned long &nbCalls);
    // true unitigs version
    void heuristic_most_covered_path_unitigs(Direction dir, Node& startingNode, Node& endNode, 
                                    int traversal_depth, std::set<Node>& usedNode, int& success, std::vector<int>& unitigs_lengths, std::vector<int>& unitigs_abundances, double& mean_abundance,
                                    unsigned int backtrackingLimit, Node *avoidFirstNode, 
                                    bool most_covered, unsigned long &nbCalls);



    std::vector<bool> interestingNodes;
};

/********************************************************************************/

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif 

