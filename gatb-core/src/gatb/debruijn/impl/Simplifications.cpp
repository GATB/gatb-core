/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014-2015  INRIA
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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_SIMPLIFCPP_
#define _GATB_CORE_DEBRUIJN_IMPL_SIMPLIFCPP_


/********************************************************************************/
// We include required definitions
/********************************************************************************/

#define DEBUG(a)    //a
#define DEBUG_TIPS(a)    //a
#define DEBUG_BULGES(a)    //a
#define DEBUG_EC(a)   //a

// the only time when you don't want to define this, is when debugging with gdb, because can't debug lambda's
#define SIMPLIFICATION_LAMBDAS 

// this is to control whether we instrument code for timing or not (shouldn't affect performance, in principle)
#define TIME(a)   a

#include <cassert>
#include <gatb/debruijn/impl/Simplifications.hpp>
#include <gatb/debruijn/impl/NodesDeleter.hpp>
#include <gatb/tools/misc/impl/Progress.hpp> // for ProgressTimerAndSystem

#include <chrono>
#define get_wtime() chrono::system_clock::now()
#define diff_wtime(x,y) (unsigned long)chrono::duration_cast<chrono::nanoseconds>(y - x).count()

#define DIR2STR(dir) ((dir==DIR_OUTCOMING) ? "outcoming" : "incoming")

using namespace std;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::system; // for System
using namespace gatb::core::system::impl; 
using namespace gatb::core::tools::dp; // for dispatcher
using namespace gatb::core::tools::dp::impl;



/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/


static const char* simplprogressFormat0 = "removing tips,    pass %2d ";
static const char* simplprogressFormat2 = "removing bulges,  pass %2d ";
static const char* simplprogressFormat3 = "removing ec,      pass %2d ";


template<typename GraphType, typename Node, typename Edge>
Simplifications<GraphType, Node, Edge>::Simplifications(GraphType *graph, int nbCores, bool verbose)
        : _nbTipRemovalPasses(0), _nbBubbleRemovalPasses(0), _nbBulgeRemovalPasses(0), _nbECRemovalPasses(0), _graph(*graph), 
        _nbCores(nbCores), _firstNodeIteration(true), _verbose(verbose)
{
    // by default; do everything
    _doTipRemoval = _doBulgeRemoval = _doECRemoval = true;

    if (graph) // may be called with graph==null in order just to get parameters
    {
        // this is just to get number of nodes
        ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem> itNode (this->_graph.iterator(), "");
        nbNodes = itNode.size();
    }

    // compute a fair amount of tips/bubble/ec after which it's useless to do another pass
    // (before, the previous system was to do a fixed amount of passes)

    cutoffEvents = std::max((uint64_t)(nbNodes / 10000), (uint64_t)1); 
    // maybe todo: estimate better and illustrate on sample genomes
    
    
    // set some default parameters for aggressive graph simplifications (minia)

    //  tips
    _tipLen_Topo_kMult = 2.5;
    _tipLen_RCTC_kMult = 10;
    _tipRCTCcutoff = 2; // SPAdes-like
    
    // bulges
    _bulgeLen_kMult = 3; 
    _bulgeLen_kAdd = 100;
    _bulgeAltPath_kAdd = 50;
    _bulgeAltPath_covMult = 1.1;

    // EC
    _ecLen_kMult = 9;
    _ecRCTCcutoff = 4;
}


/* this functions performs many rounds of all available graph simplifications 
 * this is what Minia does by default */
template<typename GraphType, typename Node, typename Edge>
void Simplifications<GraphType,Node,Edge>::simplify()
{
    unsigned long nbTipsRemoved = 0, nbTipsRemovedPreviously = 0;
    unsigned long nbBubblesRemoved = 0, nbBubblesRemovedPreviously = 0;
    unsigned long nbECRemoved = 0, nbECRemovedPreviously = 0;

    tipRemoval = "";
    bubbleRemoval = "";
    ECRemoval = "";
    
    if (_doTipRemoval)
    {
        do
        {
            nbTipsRemovedPreviously = nbTipsRemoved;
            nbTipsRemoved = removeTips();
            if (tipRemoval.size() != 0)
                tipRemoval += " + ";
            tipRemoval += to_string(nbTipsRemoved);
        }
        while ( ((nbTipsRemovedPreviously == 0 && nbTipsRemoved > 0) || (_nbTipRemovalPasses <= 2 || nbTipsRemoved >= cutoffEvents)) 
                && _nbTipRemovalPasses < 20);
    }

    if (_doBulgeRemoval)
    {
        do
        {
            nbBubblesRemovedPreviously = nbBubblesRemoved;
            nbBubblesRemoved = removeBulges(); // now we're using bulges removal, not bubbles (to follow SPAdes)
            if (bubbleRemoval.size() != 0)
                bubbleRemoval += " + ";
            bubbleRemoval += to_string(nbBubblesRemoved);
        }
        while (((nbBubblesRemovedPreviously == 0 && nbBubblesRemoved > 0) || (_nbBulgeRemovalPasses <= 2 || nbBubblesRemoved >= cutoffEvents))
                && _nbBulgeRemovalPasses < 20);
    }

    if (_doECRemoval)
    {
        do
        {
            nbECRemovedPreviously = nbECRemoved;
            nbECRemoved = removeErroneousConnections(); // now we're using bulges removal, not bubbles (to follow SPAdes)
            if (ECRemoval.size() != 0)
                ECRemoval += " + ";
            ECRemoval += to_string(nbECRemoved);
        }
        while (((nbECRemovedPreviously == 0 && nbECRemoved > 0 ) || (_nbECRemovalPasses <= 2 || nbECRemoved >= cutoffEvents))
                && _nbECRemovalPasses < 20);
    }
    
    // final simplifications rounds: do a mix of everything
    if (_doTipRemoval && _doBulgeRemoval && _doECRemoval)
    {
        nbECRemoved = 0; // reset EC removal counter
        do
        {
            nbTipsRemoved = removeTips();

            nbBubblesRemovedPreviously = nbBubblesRemoved;
            nbBubblesRemoved = removeBulges();

            nbECRemovedPreviously = nbECRemoved;
            nbECRemoved = removeErroneousConnections();

            tipRemoval += " + " + to_string(nbTipsRemoved);

            bubbleRemoval += " + " + to_string(nbBubblesRemoved);

            ECRemoval += " + " + to_string(nbECRemoved);

        }
        while (((nbECRemovedPreviously == 0 && nbECRemoved > 0) || (nbECRemoved >= cutoffEvents || nbTipsRemoved >= cutoffEvents || nbBubblesRemoved >= cutoffEvents))
                && _nbTipRemovalPasses < 30);
    }
}


// gets the mean abundance of neighboring paths around a branching node (excluding the path that starts with nodeToExclude, e.g. the tip itself)
// only considers the first 100 kmers of neighboring paths
//
template<typename GraphType, typename Node, typename Edge>
double Simplifications<GraphType,Node,Edge>::getMeanAbundanceOfNeighbors(Node& branchingNode, Node nodeToExclude)
{
    GraphVector<Edge> neighbors = _graph.neighborsEdge(branchingNode);
    unsigned int nbNeighbors = 0;
    double meanNeighborsCoverage = 0;
    //DEBUG(cout << endl << "called getMeanAbudanceOfNeighbors for node " << _graph.toString(branchingNode) << " of degrees " << _graph.indegree(branchingNode) <<"/"<< _graph.outdegree(branchingNode)<< " excluding node  " <<  _graph.toString (nodeToExclude) << endl);
    for (size_t i = 0; i < neighbors.size(); i++)
    {
        double simplePathCoverage = 0;
        if (_graph.simplePathLength(neighbors[i].from, neighbors[i].direction) > 0) // special case, we're on a unitig, don't bother querying a middle unitig kmer 
        {
            simplePathCoverage = _graph.simplePathMeanAbundance(neighbors[i].from, neighbors[i].direction);
        }
        else
        {
            Node neighbor = neighbors[i].to;
            if (neighbor == nodeToExclude) 
            {
                //DEBUG(cout << endl << "good, seen the node to exclude" << endl);
                continue; 
            }
            //std::cout << "branchingNode " << _graph.toString(neighbors[i].from) << " neighbor "  << _graph.toString(neighbor) << std::endl;;
            simplePathCoverage = _graph.simplePathMeanAbundance(neighbor, neighbors[i].direction);
        }
        meanNeighborsCoverage += simplePathCoverage;
        nbNeighbors++;
        //DEBUG(cout << endl << "got simple path coverage for neighbor " << nbNeighbors  << " : " << " meancoverage: " <<simplePathCoverage << " over " << pathLen << " kmers" << endl);
    }
    meanNeighborsCoverage /= nbNeighbors;
    return meanNeighborsCoverage;
}

// this needs to be in Graph.cpp of gatb-core
template<typename GraphType, typename Node, typename Edge>
string Simplifications<GraphType,Node,Edge>::path2string(Direction dir, Path_t<Node> p, Node endNode)
{
    // naive conversion from path to string
    string p_str;
    if (dir == DIR_INCOMING)
    {
        p_str = "";
        for (size_t i = 0; i < p.size(); i++)
            p_str.push_back(p.ascii(p.size()-1-i));
        p_str += _graph.toString(p.start);
    }
    else
    {
        p_str = _graph.toString(p.start);
        for (size_t i = 0; i < p.size(); i++)
            p_str.push_back(p.ascii(i));
    }
    //std::cout << "path2string " << to_string(dir) << " " << p_str <<endl;
    return p_str;
}

// in the case of the new heuristic_most_covered, we only record the nucleotide at branching nodes
// this function will only be used for debugging the paths
template<typename Node>
string unitig2string(Direction dir, Path_t<Node> p, Node endNode)
{
    return "";
}
 

// this needs to be in Graph.cpp of gatb-core
// computes mean abundance of path
// optional params!
// skip_first: skip N nucleotides at beginning of path
// skip_last: skip N nucleotides at end of path
// e.g. skip_first=skip_last=1 will skip the first and last kmer of path
template<typename GraphType, typename Node, typename Edge>
double Simplifications<GraphType,Node,Edge>::path2abundance(Direction dir, Path_t<Node> p, Node endNode, unsigned int skip_first, unsigned int skip_last)
{
    // the buildnode call here needs to be adapted for NodeGU
#if 0
    string s = path2string(dir,p,endNode);
    if (p.size() == 0) return 0;
    if (p.size() <= skip_first + skip_last) return 0;

    unsigned long mean_abundance = 0;
    //std::cout << " abundance: " ;
    for (size_t i = skip_first; i < p.size() - skip_last; i++)
    {
        Node node = _graph.buildNode((char *)(s.c_str()), i);
        int abundance = _graph.queryAbundance(node);
        mean_abundance += abundance;
        //std::cout << " " << to_string(abundance);
        //cout << endl << "node: " << _graph.toString (node) << " abundance: "  << to_string(abundance) << endl;
    }                              
    //std::cout << std::endl;
    mean_abundance /= (p.size() - skip_first - skip_last);
    return mean_abundance;
#endif
    return 0;
}

inline string maybe_print(long value, string str)
{
    if (value == 0)
        return "";
    return to_string(value) + " " + str;
}


/* the RCTC terms is taken from SPADES and means Relative Coverage Tip Clipping
 * but here i'm using it loosely to mean "Relative Coverage"
 * this function gets coverage of the simple path,
   then compares it to coverage of other paths connected to the last node of it. */
template<typename GraphType, typename Node, typename Edge>
bool Simplifications<GraphType,Node,Edge>::satisfyRCTC(double pathAbundance, Node& node, double RCTCcutoff, Direction dir)
{
    // explore the other two or more simple paths connected to that path, to get an abundance estimate
    // but first, get the branching node(s) the tip is connected to 
    //
    // invariant: node is the last node of a path. 

    GraphVector<Edge> connectedBranchingNodes = _graph.neighborsEdge(node, dir);
    unsigned int nbBranchingNodes = 0;
    double meanNeighborsCoverage = 0;

    // disabling the following check, becuse in EC removal, reverse direction, node may not be a last node (was in a hurry, need to revisit whole RCTC someday anyway to see if it can be improved)
/*    if (_graph.simplePathLength(node , dir) > 0) 
    {
        std::cout << "satisfyRCTC; unexpected: node isn't a last node?" << std::endl;
        std::cout << _graph.toString(node) << " dir " << dir << " simple path length: " << _graph.simplePathLength(node , dir)<< " connectedBranchingNodes[0].to: " << _graph.toString(connectedBranchingNodes[0].to) << " connectedBranchingNodes.size(): " << connectedBranchingNodes.size() << std::endl;
        exit(1);
    }*/

    /* there may be more than one connectedBranchNode
     * it's weird when it's more than one branching node though, it's a situation like:
     * 
     *                ..o--o--o--o..
     *                    /
     *  nodes-> o--o--o--o
     *                    \
     *                ..o--o--o--o.. 
     *
     *  instead of the classical:
     *
     *  nodes-> o--o--o--o
     *                    \
     *                ..o--o--o--o..)*/

    for (size_t j = 0; j < connectedBranchingNodes.size(); j++)
    {
        meanNeighborsCoverage += this->getMeanAbundanceOfNeighbors(connectedBranchingNodes[j].to, node);
        nbBranchingNodes++;
    }
    if (nbBranchingNodes > 0)
        meanNeighborsCoverage /= nbBranchingNodes;

    bool isRCTC = (meanNeighborsCoverage > RCTCcutoff * pathAbundance);

    //DEBUG(cout << endl << "RCTC test, over " << nbBranchingNodes << " connected nodes. Global mean neighbors coverage: " << meanNeighborsCoverage <<  " compared to mean tip abundance over "<< nodes.size() << " values : " << meanTipAbundance << (debugstdev? " stddev: " : "") << (debugstdev? to_string(stdevTipAbundance): "") << ", is RCTC satisfied? " << isRCTC << endl);

    return isRCTC;
}

/* this function removes tips in the graph, using an algorithm designed in SPAdes
 *
 * here is an in-depth analysis of SPAdes 3.5 tip clipping conditions: (following graph_simplifications.hpp and simplifications.info)
 *
 * * tc_lb is a coefficient, setting tip length to be max(read length, g*tc_lb) (see simplification_settings.hpp);
 *
 * * cb is a plain coverage upper bound (see basic_edge_conditions.hpp)
 * when it's auto, it's equal to ec_bound. (source: set_detected_coverage_bound(gp.ginfo.ec_bound() in graph_simplification.hpp)
 *
 * * * ec_bound is given by the error threshold found by the coverage model in kmer_coverage_model.cpp (see: genomic_info_filler.cpp)
 * it's the largest abundance value such that the probability of an erroneous kmer is less than 0.05. And, if it's smaller than the valley, it's the valley.
 * EC bound looks off in metagenomic assemblies (a value over 50.. seems high).
 *
 * * * in single cell mode, it's max(average edge coverage, graph theshold), and graph threshold is quite advanced! (omni_tools.hpp)
 * but already, in multicell mode, spades does a good job at tip clipping, so let's look at that first.
 *
 * * rctc is a "relative coverage tip condition" and the value given next to it is: max_relative_coverage
 * rctc examines whether the coverage of the tip is <= max_relative_coverage * (max_{over all neighbors}(coverage of neighbor) + 1)
 
 * actual tip clipping code for "tc", i.e. implementations of the conditions, is in tip_clipper
 * (omni/graph_processing_algorithm.hpp:EdgeRemovingAlgorithm() is just generically parsing conditions)
 *
 * in one condition, tc_lb = 3.5, cb is 1000000 (minia's coverages values don't go that far.. yet..), rctc 2
 * and in another, tc_lb is 10 and cb is auto
 *
 * some historical facts:
 * for CAMI we were more loose than SPAdes: topological tip had no cov criterion, wherehas it should have had rctc (like spades)
 * and long tc_lb10 tips should have the auto coverage bound, but instead they had rctc 2. 
 *
 * here we also keep that philosophy: topological tips (<=3.5*k) do not need a coverage criterium to be removed (no rctc)
 *
 * so TODO: make it more strict. but for now I'm focusing on EC.
 */
template<typename GraphType, typename Node, typename Edge>
unsigned long Simplifications<GraphType,Node,Edge>::removeTips()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int _maxTipLengthTopological = (unsigned int)((float)k * _tipLen_Topo_kMult); // aggressive with SPAdes length threshold, but no coverage criterion
    unsigned int _maxTipLengthRCTC = (unsigned int)(k * _tipLen_RCTC_kMult); // experimental, SPAdes-like

    unsigned long nbTipsRemoved = 0;

    // stats
    unsigned long timeAll = 0, timeDecision = 0, timeProcessing = 0, timeSimplePath = 0, timeSimplePathLong = 0, timeSimplePathShortTopo = 0, timeSimplePathShortRCTC = 0, timeDel = 0, timeCache = 0;
    unsigned long nbTipCandidates = 0;


    char buffer[128];
    sprintf(buffer, simplprogressFormat0, ++_nbTipRemovalPasses);
    /** We get an iterator over all nodes */
    /* in case of pass > 1, only over cached branching nodes */
    // because in later iterations, we have cached non-simple nodes, so iterate on them
    ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem> *itNode; 
    if (_firstNodeIteration )
    {
        itNode = new ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem>(_graph.GraphType::iterator(), buffer, _verbose);
        if (_verbose)
            std::cout << "iterating on " << itNode->size() << " nodes on disk" << std::endl;
    }
    else
    {
        itNode = new ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem>(_graph.GraphType::iteratorCachedNodes(), buffer, _verbose);
        if (_verbose)
            std::cout << "iterating on " << itNode->size() << " cached nodes" << std::endl;
    }

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);

    // nodes deleter stuff
    NodesDeleter<Node,Edge,GraphType> nodesDeleter(_graph, nbNodes, _nbCores, _verbose);

    dispatcher.iterate (*itNode, [&] (Node& node)
    {
         /* just a quick note, which was observed in the context of flagging some node as uninteresting (not used anymore).
          * property: "a tip (detected at some point after some rounds of simplifications) is not necessarily a branching node initially in the original graph"
         * here a strange dbg motif:
         *
         *  
         *       (t)
         *       / \
         *      /   \
         *     v     v
         *  -->o    (n1)-->
         *
         *  node n1 is simple, but node (t) is considered a tip so will be deleted.
         *  thus n1 will become a tip too
         *  question is, should we actually delete (t)?
         *  i'd argue so, it's likely artefactual. (n1) could be artefactual too.
         *  consider those erroneous kmers:
         *
         *  AAT
         *   ATC
         *   ATG
         *
         * and assuming a true genome is TCTGTC (circular)
         *
         * so graph is:
         *
         *                   /-------------\
         *                  v              |
         *   AAT -> ATC -> TCT -> CTG      |
         *    \                   /        |
         *     \                 /         |
         *      \               v          |
         *       \-> ATG  --> TGT --> GTC -/
         *
         *  the three leftmost nodes should be removed and the four rightmost nodes are correct.
         *  once AAT is removed, ATC and ATG are legit tips; but they were simple nodes initially. QED
         *
         */

        TIME(auto start_thread_t=get_wtime());

        // skip deleted nodes
        if (_graph.isNodeDeleted(node)) {
                TIME(auto end_thread_t=get_wtime()); 
                TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t))); 
            return; }  

        unsigned inDegree = _graph.indegree(node), outDegree = _graph.outdegree(node);

        /* tips have out/in degree of 0 on one side, and any non-zero degree on the other */
        if ((inDegree == 0 || outDegree == 0) && (inDegree != 0 || outDegree != 0))
        {
            bool isShortTopological = true;
            bool isShortRCTC = true;

            DEBUG_TIPS(cout << endl << "deadend node: " << _graph.toString (node) << endl);
            __sync_fetch_and_add(&nbTipCandidates,1);

            // this call is only to get the direction. a bit hacky.
            GraphVector<Edge> neighbors = _graph.neighborsEdge(node); 
            // but in fact, node may have one or more neighbors in that direction
           
            if (neighbors.size() == 0) { std::cout << "unexpected problem during removeTips, no neighbor; " << inDegree << " " << outDegree << " " << _graph.toString(node) << std::endl; exit(1);}

            /* it may appear that we're only going to follow its first neighbor, but in fact, neighbors[0].from is node */
            /* so, follow the simple path from this start tip node to the further node that has out-branching (out is w.r.t to the direction) */
            TIME(auto start_simplepath_t=get_wtime());
            Node&     simplePathStart = neighbors[0].from;
            Direction simplePathDir   = neighbors[0].direction;
            unsigned int pathLen = _graph.simplePathLength(simplePathStart,simplePathDir);
            TIME(auto end_simplepath_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeSimplePath, diff_wtime(start_simplepath_t,end_simplepath_t)));
           
            if (k + pathLen > _maxTipLengthTopological) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
               isShortTopological = false;
            if (k + pathLen > _maxTipLengthRCTC) 
               isShortRCTC = false;

            if (isShortTopological)
                TIME(__sync_fetch_and_add(&timeSimplePathShortTopo, diff_wtime(start_simplepath_t,end_simplepath_t)));
            if (isShortRCTC)
                TIME(__sync_fetch_and_add(&timeSimplePathShortRCTC, diff_wtime(start_simplepath_t,end_simplepath_t)));
    
            // if it's not short, then no point in computing whether it's connected
            // also, it's pointless to examine it later, as it will never become short again (we only delete nodes)
            // so mark the origin as non interesting! (big speed up)
            if ( ! (isShortTopological || isShortRCTC) )
            {   
                TIME(__sync_fetch_and_add(&timeSimplePathLong, diff_wtime(start_simplepath_t,end_simplepath_t)));
                TIME(auto end_thread_t=get_wtime()); 
                TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t)));
                return;
            }
            
            TIME(auto start_tip_decision_t=get_wtime());

            // at this point, the last node in "nodes" is the last node of the tip.
            // check if it's connected to something. 
            // condition: degree > 1, because connected to the tip and to that "something"
            Node lastNode           = _graph.simplePathLastNode(simplePathStart,simplePathDir);
            bool isConnected = (_graph.neighborsEdge(lastNode, simplePathDir).size() >= 1);
            if (pathLen == 0)
            {
                // special case: only a single tip node, check if it's not isolated
                isConnected |=  (_graph.indegree(node) != 0 || _graph.outdegree(node) != 0); 
            }

            // TODO would be worth it to check if the node where the tip is connected to, is also connected to the tip and something else. i.e existence of a V-shaped pattern.

            bool isTopologicalShortTip = isShortTopological && isConnected; 
            bool isMaybeRCTCTip = isShortRCTC && isConnected;

            DEBUG_TIPS(cout << endl << "pathlen: " << pathLen << " last node " << _graph.toString(lastNode) << " neighbors in/out: " <<_graph.indegree(lastNode) << " " << _graph.outdegree(lastNode) << " istoposhorttip: " << isTopologicalShortTip << endl);
            
            double pathMeanAbundance = _graph.simplePathMeanAbundance(simplePathStart,simplePathDir);

            bool isRCTCTip = false;
            if (!isTopologicalShortTip && isMaybeRCTCTip)
            {
                isRCTCTip = this->satisfyRCTC(pathMeanAbundance, lastNode, _tipRCTCcutoff, simplePathDir); 
                /* fun fact: not putting "this->" crashes gcc 4.7; was fun to debug :\ */
            }

            bool isTip = isTopologicalShortTip || isRCTCTip; 
            
            TIME(auto end_tip_decision_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeDecision, diff_wtime(start_tip_decision_t,end_tip_decision_t)));

            TIME(auto start_tip_processing_t=get_wtime());

            if (isTip)
            {
                if (nodesDeleter.get(simplePathStart))
                    {
                        // not double-counting that delete
                    }
                else
                {
                    // delete it
                    _graph.simplePathDelete(simplePathStart, simplePathDir, nodesDeleter);

                    __sync_fetch_and_add(&nbTipsRemoved, 1);
                    DEBUG_TIPS(cout << endl << "TIP FOUND, deleting node : " << _graph.toString(simplePathStart) << endl);
                }
            } // end if isTip

            TIME(auto end_tip_processing_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeProcessing, diff_wtime(start_tip_processing_t,end_tip_processing_t)));

        } // end if degree correspond to putative end-of-tip

        TIME(auto end_thread_t=get_wtime()); 
        TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t)));

    }); // parallel

    TIME(auto start_nodesdel_t=get_wtime());
    DEBUG_TIPS(std::cout << "end of tip removal pass" << std::endl;);
    
    // now delete all nodes, in parallel
    nodesDeleter.flush();

    TIME(auto end_nodesdel_t=get_wtime()); 
    TIME(__sync_fetch_and_add(&timeDel, diff_wtime(start_nodesdel_t,end_nodesdel_t)));


    TIME(auto start_nodescache_t=get_wtime());

    if (_firstNodeIteration)
        _graph.cacheNonSimpleNodes(_nbCores, true);

    TIME(auto end_nodescache_t=get_wtime()); 
    TIME(__sync_fetch_and_add(&timeCache, diff_wtime(start_nodescache_t,end_nodescache_t)));
 
    // stats
    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);
    if (_verbose)
    {
        cout << nbTipsRemoved << " tips removed. " << endl;
        cout << nbTipCandidates << " tip candidates passed degree check. " << endl;
        TIME(cout << "Tips timings:   " << timeAll / unit << " CPUsecs total."<< endl);
        TIME(cout << "                " << timeSimplePath / unit << " CPUsecs simple path traversal, including:" << endl);
        TIME(cout << "                    " << timeSimplePathLong  / unit << " CPUsecs long simple paths" << endl);
        TIME(cout << "                    " << timeSimplePathShortTopo / unit << " CPUsecs short (topological) simple paths" << endl);
        TIME(cout << "                    " << timeSimplePathShortRCTC / unit << " CPUsecs short (RCTC) simple paths" << endl);
        TIME(cout << "                " << timeDecision   / unit << " CPUsecs tip decision" << endl);
        TIME(cout << "                " << timeProcessing / unit << " CPUsecs tip processing" << endl);
        TIME(cout << "Nodes deletion:   " << timeDel / unit << " CPUsecs."<< endl);
        TIME(cout << "Nodes caching :   " << timeCache / unit << " CPUsecs."<< endl);
    }
    
    _firstNodeIteration = false;

    return nbTipsRemoved;
}

enum HMCP_Success { HMCP_DIDNT_FIND_END = 0, HMCP_FOUND_END = 1 , HMCP_MAX_DEPTH = -1, HMCP_LOOP = -2};

static string hmcpstatus2ascii(int h){
    switch (h)
    {
        case HMCP_DIDNT_FIND_END:
            return "didn't find end";
        case HMCP_FOUND_END:
            return "success, found end";
        case HMCP_MAX_DEPTH:
            return "max depth reached";
        case HMCP_LOOP:
            return "loop found, aborted";
    }
    return "unknown";
}


static double unitigs_chain2abundance(vector<int> &unitigs_lengths, vector<int> &unitigs_abundances)
{
    assert(unitigs_lengths.size() == unitigs_abundances.size());
    //std::cout << "unitigs lengths " << unitigs_lengths.size() << std::endl;
    double mean_abundance = 0;
    int sum_lengths = 0;
    for (unsigned int i = 0; i < unitigs_lengths.size(); i++)
    {
        mean_abundance += unitigs_lengths[i] * unitigs_abundances[i];
        sum_lengths += unitigs_lengths[i];
        //std::cout << "mean abundance now: "  << to_string(mean_abundance) << std::endl;
        //std::cout << "sum lengths now: "  << to_string(sum_lengths) << std::endl;
    }
    return mean_abundance / sum_lengths;
}


/* this function finds the most covered path between start and end node
 * note: the returned mean abundance does not include start and end nodes */
// endNode is a node just after the bulge path. it's branching.
template<typename GraphType, typename Node, typename Edge>
void Simplifications<GraphType,Node,Edge>::heuristic_most_covered_path(
        Direction dir, Node& startNode, Node& endNode, 
        int traversal_depth, int& success, double &mean_abundance, 
        Path_t<Node> &res_path,
        unsigned int backtrackingLimit, Node *avoidFirstNode,
        bool most_covered, bool kmer_version /* was the original version, where nodes were traversed kmer-by-kmer, now it's one simple path at a time */)
{
    set<Node> usedNode;
    usedNode.insert(startNode);
    Path_t<Node> current_path;
    current_path.start = startNode;
    success = HMCP_DIDNT_FIND_END;
    vector<int> abundances; 
    vector<int> unitigs_lengths, unitigs_abundances; 
    unsigned long nbCalls = 0;
    mean_abundance = 0;

    // so i've disabled the kmer version, see comment on hmcp_old() function below
#if 0
    if (kmer_version)
    {
        heuristic_most_covered_path_old(dir, startNode, endNode, traversal_depth, current_path, usedNode, success, abundances,
                backtrackingLimit, avoidFirstNode, 
                most_covered, 
                res_path,
                nbCalls);

        if (success == HMCP_FOUND_END)
        { // no need to compute abundance in failed cases
            for (unsigned int i = 0; i < abundances.size(); i++){
                mean_abundance += abundances[i];
            }
            mean_abundance /= abundances.size();
        }
    }
    else
    {
        // intermediate version, before GraphUnitigs
/*        heuristic_most_covered_path(dir, startNode, endNode, traversal_depth, current_path, usedNode, success, mean_abundance,
                backtrackingLimit, avoidFirstNode, 
                most_covered, 
                res_path,
                nbCalls);
  */          
        
        // WARNING:this is slow if ran on Graph. but fast on GraphUnitigs.
#endif    
        heuristic_most_covered_path_unitigs(dir, startNode, endNode, 
                traversal_depth, usedNode, success, unitigs_lengths, unitigs_abundances, mean_abundance,
                backtrackingLimit, avoidFirstNode, 
                most_covered, nbCalls);
#if 0
    }
#endif

    /* below this point: debug stuff*/

    bool debug_abundances = false;
    if (debug_abundances)
    {
        cout << "abundance for path (is most covered path: " << most_covered << "): ";
        for (unsigned int i = 0; i < abundances.size(); i++)
            cout << abundances[i]<< " ";
        cout << ";"<<endl;
    }

    bool debug_nbcalls = false;
    if (debug_nbcalls)
        cout << "number of path-finding calls: " << nbCalls << endl;
}

#if 0 //the nt stuff is not supported in EdgeGU anymore, so i'm disabling this code

template<typename GraphType, typename Node, typename Edge>
void Simplifications<GraphType,Node,Edge>::heuristic_most_covered_path_old(
        Direction dir, Node& startNode, Node& endNode,
        int traversal_depth, Path_t<Node>& current_path, set<Node>& usedNode, int& success, vector<int>& abundances,
        unsigned int backtrackingLimit, Node *avoidFirstNode, 
        bool most_covered, Path_t<Node> &res_path,
        unsigned long &nbCalls)
{
    // inspired by all_consensuses_between
    nbCalls++;
    
    if (traversal_depth < -1)
    {
        success = HMCP_DIDNT_FIND_END;
        return;
    }

    if (startNode == endNode)
    {
        success = HMCP_FOUND_END;
        res_path = current_path;
        return;
    }

    GraphVector<Edge> neighbors = _graph.neighborsEdge (startNode, dir);

    /** We loop these neighbors. */
    vector<std::pair<int, Edge> > abundance_node;
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& edge = neighbors[i];

        if (avoidFirstNode != NULL && edge.to == *avoidFirstNode)
            continue;

        // don't resolve bubbles containing loops
        // (tandem repeats make things more complicated)
        // that's a job for a gapfiller
        if (usedNode.find(edge.to) != usedNode.end())
        {
            success = HMCP_LOOP;
            return;
        }
        
        unsigned int abundance = _graph.queryAbundance(neighbors[i].to);
        abundance_node.push_back(std::make_pair(abundance, edge));
        //cout << endl << "node: " << _graph.toString (neighbors[i].to) << " abundance: "  << to_string(abundance) << endl;
    }

    std::sort(abundance_node.begin(), abundance_node.end()); // sort nodes by abundance
    if (most_covered) 
        std::reverse(abundance_node.begin(), abundance_node.end()); // decreasing abundances

    // traverse graph in abundance order, return most abundant path
    for (unsigned int i = 0; i < abundance_node.size(); i++)
    {
        Edge edge = abundance_node[i].second;

        // generate extended consensus sequence
        Path_t<Node> extended_path(current_path);
        extended_path.push_back (edge.nt);

        // generate list of used kmers (to prevent loops)
        set<Node> extended_kmers (usedNode);
        extended_kmers.insert (edge.to);

        // extend abundances
        vector<int> extended_abundances (abundances);
        if (edge.to != endNode) // skip abundance of last node
            extended_abundances.push_back(abundance_node[i].first);

        // recursive call to all_consensuses_between
        heuristic_most_covered_path_old (
            dir,
            edge.to,
            endNode,
            traversal_depth - 1,
            extended_path,
            extended_kmers,
            success,
            extended_abundances,
            backtrackingLimit,
            NULL, // no longer avoid nodes
            most_covered,
            res_path,
            nbCalls
        );

        if (backtrackingLimit > 0 && nbCalls >= backtrackingLimit)
        {
            success = HMCP_MAX_DEPTH;
            return;
        }

        if (success == HMCP_FOUND_END)
        {
            abundances = extended_abundances;
            return; 
        }
    }
    return;
}

/* faster variant of the algo above, with cached simple paths instead of traversing node-by-node */
// Note: not totally equivalent to the kmer version, some bubbles are found in the kmer version but not in this one
template<typename GraphType, typename Node, typename Edge>
void Simplifications<GraphType,Node,Edge>::heuristic_most_covered_path(
        Direction dir, Node& startNode, Node& endNode, 
        int traversal_depth, Path_t<Node>& current_path, set<Node>& usedNode, int& success, double &mean_abundance,
        unsigned int backtrackingLimit, Node *avoidFirstNode, 
        bool most_covered, Path_t<Node> &res_path,
        unsigned long &nbCalls)
{
    nbCalls++;
    
    if (traversal_depth < -1)
    {
        success = HMCP_DIDNT_FIND_END;
        return;
    }

    Node current_node = startNode;
 
    if (current_node == endNode)
    {
        success = HMCP_FOUND_END;
        res_path = current_path;
        mean_abundance = path2abundance(dir, res_path, endNode, 1, 1);
        // to get the seq, make some sort of unitig2string(), but i never did it
        return;
    }

    set<Node>& traversedNodes (usedNode);

    Path_t<Node> current_extended_path(current_path);
    int extra_depth = 1;

    auto extend_and_check = [&](Edge &edge)
    {
        current_node = edge.to;
        current_extended_path.push_back( edge.nt );
        extra_depth++;
        nbCalls++;

        if (traversedNodes.find(current_node) != traversedNodes.end() )  // loop
        {
            success = HMCP_LOOP;
            return true;
        }
 
        traversedNodes.insert(current_node); 

        if (current_node == endNode)
        {
            success = HMCP_FOUND_END;
            res_path = current_extended_path;
            mean_abundance = path2abundance(dir, res_path, endNode, 1, 1);
            return true;
        }

        return false;
    };


    // traverse simple path from that node
    // we end up at a branching node. 
    // if the node has no out-branching, keep going, we don't care about in-branching here.
    do
    {
        GraphIterator <Edge> itEdges
            = _graph.template simplePathEdge (current_node, dir);
        
        for (itEdges.first(); !itEdges.isDone(); itEdges.next())
        {
            if (extend_and_check(*itEdges)) // updates current_node and current_extended_path, as well as nbCalls and extra_depth
                return;
        }

        // end of simple path, yet no out-branching? means there is in-branching
        if (_graph.degree(current_node, dir) == 1)
        {
            // get the node after in-branching
            GraphVector<Edge> neighbors = _graph.neighborsEdge(current_node, dir);
        
            if (extend_and_check(neighbors[0])) // updates current_node and current_extended_path, as well as nbCalls and extra_depth
                return;
        }
        else
            break; // we're either at a deadend or at a branching node
    }
    while (true);
            
    // at this point, we've traversed simple paths from startNode and we ended up at a branching node or a deadend

    // get neighbors of branching node
    GraphVector<Edge> neighbors = _graph.neighborsEdge (current_node, dir);

    /** We loop over the neighbors of that branching node, to order them by abundance */
    vector<std::pair<int, Edge> > abundance_node;
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& edge = neighbors[i];

        if (avoidFirstNode != NULL && edge.to== *avoidFirstNode)
            continue;

        // don't resolve bubbles containing loops
        // (tandem repeats make things more complicated)
        // that's a job for a gapfiller
        if (traversedNodes.find(edge.to) != traversedNodes.end())
        {
            success = HMCP_LOOP;
            return;
        }
        
        unsigned int abundance = _graph.queryAbundance(neighbors[i].to);
        abundance_node.push_back(std::make_pair(abundance, edge));
    }

    std::sort(abundance_node.begin(), abundance_node.end()); // sort nodes by abundance
    if (most_covered) 
        std::reverse(abundance_node.begin(), abundance_node.end()); // decreasing abundances

    // traverse graph in abundance order, return most abundant path
    for (unsigned int i = 0; i < abundance_node.size(); i++)
    {
        Edge edge = abundance_node[i].second;

        // generate extended consensus sequence
        Path_t<Node> extended_path(current_extended_path);
        extended_path.push_back (edge.nt);

        // generate list of used kmers (to prevent loops)
        set<Node> extended_kmers (traversedNodes);
        extended_kmers.insert (edge.to);

        // recursive call to all_consensuses_between
        heuristic_most_covered_path (
            dir,
            edge.to,
            endNode,
            traversal_depth - extra_depth,
            extended_path,
            extended_kmers,
            success,
            mean_abundance,
            backtrackingLimit,
            NULL, // no longer avoid nodes
            most_covered,
            res_path,
            nbCalls
        );
              
        if (backtrackingLimit > 0 && nbCalls >= backtrackingLimit)// if no more backtracking, return immediately
        {
            success = HMCP_MAX_DEPTH;
            return;
        }

        if (success == HMCP_FOUND_END)
        {
            return; 
        }
 
    }

    return;
}
#endif

/* GraphUnitigs version of the algo above. 
 * doesn't construct a Path anymore, and handles abundances differently, so i preferred to make another function rather than re-using the one above */
template<typename GraphType, typename Node, typename Edge>
void Simplifications<GraphType,Node,Edge>::heuristic_most_covered_path_unitigs(
        Direction dir, Node& startNode, Node& endNode, 
        int traversal_depth, set<Node>& usedNode, int& success, vector<int>& unitigs_lengths, vector<int>& unitigs_abundances, double& mean_abundance,
        unsigned int backtrackingLimit, Node *avoidFirstNode, 
        bool most_covered, unsigned long &nbCalls)
{
    bool debug = false;
    
    if (traversal_depth < -1)
    {
        success = HMCP_DIDNT_FIND_END;
        return;
    }
    
    nbCalls++;

    Node current_node = startNode;
    if (debug)
        std::cout << "HMCP rec, now at : " << _graph.toString(current_node) << " dir " << dir<< std::endl;;
 
    if (current_node == endNode)
    {
        success = HMCP_FOUND_END;
        mean_abundance = unitigs_chain2abundance(unitigs_lengths, unitigs_abundances);
        return;
    }

    set<Node>& traversedNodes (usedNode);
    int extra_depth = 1;

    auto processNode = [&unitigs_lengths, &unitigs_abundances, &endNode, &current_node, &success, &mean_abundance, &traversedNodes](Node &node)
    {
        current_node = node;
        if (current_node == endNode)
        {
            success = HMCP_FOUND_END;
            mean_abundance = unitigs_chain2abundance(unitigs_lengths, unitigs_abundances);
            return true;
        }
        if (traversedNodes.find(current_node) != traversedNodes.end() )  // loop
        {
            //std::cout << "loop : " << _graph.toString(current_node) << std::endl;;
            success = HMCP_LOOP;
            return true;
        }
        traversedNodes.insert(current_node); 
        return false;
    };

    // traverse simple path from that node
    // we end up at a branching node. 
    // if the node has no out-branching, keep going, we don't care about in-branching here.
    do
    {
        if (debug)
            std::cout << "HMCP traversing simple path from node: " << _graph.toString(current_node) << std::endl;;
        
        //if (processNode(current_node)) // don't check whether we're done here, because first kmer is already inserted in traversedKmers by loop invariant

        Node&     simplePathStart = current_node;
        Direction simplePathDir   = dir;
        unsigned int pathLen = _graph.simplePathLength(simplePathStart,simplePathDir);
        if (pathLen > 0)
        { 
            double pathMeanAbundance = _graph.simplePathMeanAbundance(simplePathStart,simplePathDir);
            
            unitigs_lengths.push_back(pathLen);
            unitigs_abundances.push_back(pathMeanAbundance);

            nbCalls += pathLen + 1;
            extra_depth += pathLen + 1; 
        
            Node lastNode           = _graph.simplePathLastNode     (simplePathStart,simplePathDir);
            if (processNode(lastNode)) // verify whether we're done
                return;

            if (debug)
                std::cout << "HMCP now at last node : " << _graph.toString(current_node) << std::endl;;


        }
        else
        {
            if (debug)
            {
                Node lastNode           = _graph.simplePathLastNode     (simplePathStart,simplePathDir);
                std::cout << "HMCP last node was equal to first node: " << _graph.toString(lastNode) << " " << _graph.toString(current_node) << std::endl;
            }
        }

        // end of simple path, yet no out-branching? means there is in-branching
        if (_graph.degree(current_node, dir) == 1)
        {
            // get the node after in-branching
            GraphVector<Edge> neighbors = _graph.neighborsEdge(current_node, dir);

            current_node = neighbors[0].to;
            if (_graph.degree(current_node, reverse(dir)) <= 1)
            {
                std::cout << "Weird, there was supposed to be an in-neighbor. Maybe there's a loop. Remove this print if it never happens" << std::endl;
                return;
            }

             if (debug)
                std::cout << "end of simple path, at node, there is in-branching: " << _graph.toString(current_node) << " degrees " << _graph.degree(current_node, reverse(dir)) << "/" << _graph.degree(current_node, dir) << std::endl;;


            nbCalls ++;
            extra_depth++; 
        
            if (processNode(current_node)) // verify whether we're done
                return;
        }
        else
            break; // we're either at a deadend or at a branching node
    }
    while (true);
            
    // at this point, we've traversed simple paths from startNode and we ended up at a branching node or a deadend

    // get neighbors of branching node
    GraphVector<Edge> neighbors = _graph.neighborsEdge (current_node, dir);
    
    assert(neighbors.size() != 1);
    
    /** We loop over the neighbors of that branching node, to order them by abundance 
     * there's a variant here from other flavors of HMP: we take the most abundance path. */
    vector<std::pair<int, Edge> > abundance_node;
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& edge = neighbors[i];

        if (avoidFirstNode != NULL && edge.to == *avoidFirstNode)
            continue;

        // don't resolve bubbles containing loops
        // (tandem repeats make things more complicated)
        // that's a job for a gapfiller
        if (traversedNodes.find(edge.to) != traversedNodes.end())
        {
            success = HMCP_LOOP;
            return;
        }
        
        double abundance = _graph.simplePathMeanAbundance(edge.to,dir);
        abundance_node.push_back(std::make_pair(abundance, edge));
    }

    std::sort(abundance_node.begin(), abundance_node.end()); // sort nodes by abundance
    if (most_covered) 
        std::reverse(abundance_node.begin(), abundance_node.end()); // decreasing abundances

    // traverse graph in abundance order, return most abundant path
    for (unsigned int i = 0; i < abundance_node.size(); i++)
    {
        Edge edge = abundance_node[i].second;

        // generate list of used kmers (to prevent loops)
        set<Node> extended_kmers (traversedNodes);
        extended_kmers.insert (edge.to);

        vector<int> new_unitigs_lengths(unitigs_lengths);
        vector<int> new_unitigs_abundances(unitigs_abundances);

        // recursive call to all_consensuses_between
        heuristic_most_covered_path_unitigs(
            dir,
            edge.to,
            endNode,
            traversal_depth - extra_depth,
            extended_kmers,
            success,
            new_unitigs_lengths,
            new_unitigs_abundances,
            mean_abundance,
            backtrackingLimit,
            NULL, // no longer avoid nodes
            most_covered,
            nbCalls
        );
              
        if (backtrackingLimit > 0 && nbCalls >= backtrackingLimit)// if no more backtracking, return immediately
        {
            success = HMCP_MAX_DEPTH;
            return;
        }

        if (success == HMCP_FOUND_END)
        {
            return; 
        }
 
    }

    return;
}

class DebugBR
{
    string start;
    struct InfoNeighbor
    {
        int idx;
        string neighbor;
        string dir;
        string lastNode;
        int pathLen;
        int out_degree_lastNode;
        string branchTo;
        bool isTopologicalBulge;
    };
    vector<InfoNeighbor> neighbors;
   
    public: 
    DebugBR(string node, int nb_neighbors) : start(node), neighbors(nb_neighbors) { 
        // std::cout << "new candidate bulge: " << node << std::endl;
    }
    
    void infoNeighbor(int neighbor_idx, string neighbor, string direction, string lastNode, int pathLen, int nb_outneighbors, string branchTo, bool isTopologicalBulge)
    {   
        InfoNeighbor info{neighbor_idx, neighbor, direction, lastNode, pathLen, nb_outneighbors, branchTo, isTopologicalBulge};
        neighbors[neighbor_idx] = info;
        //cout << "last node of simple path: "<< lastNode << " and " << nb_outneighbors << " neighbors in bubble direction" << endl;
        //cout << "endNode: " << branchTo << endl);
        //cout << "pathlen: " << pathLen << " istopobulge: " << isTopologicalBulge << endl);
    }

    void draw()
    {
       cout << start << " candidate bulge" << std::endl;;
       for (unsigned int i = 0; i < neighbors.size(); i++)
       {
            const InfoNeighbor &nfo = neighbors[i];
            std::cout << " --> " << nfo.neighbor << "(" << nfo.dir << ") ...[len=" << nfo.pathLen << "]... " << nfo.lastNode << " [neighbors: " << nfo.out_degree_lastNode << "] -> " << nfo.branchTo << " istopobulge: " << nfo.isTopologicalBulge << std::endl;
       }
    }
};

/* bulge removal algorithm. mimics spades, which doesnt remove bubbles, but only bulges. looks as effective.
 *
 * see a-b-c here for an explanation of bulge removal: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791033/figure/f4/
 *
 * spades pops bulges based on something like the ratio between most examined simple path and a more covered path is (whether it is above 1.1).
 * so i'm actually doing just that. I recall checking spades source code to implement this. this was during CAMI.
 *
 * Small apparte about simple paths: in SPAdes' source, a simple path isn't a non-branching one, it is rather the wikipedia definition: one where nodes aren't repeated (and also here, no node is its own reverse-complement). This makes me think that GATB's simplePath function is a bit of a misnomer, should be called nonBranchingPath.
 *
 * so.. this is what we did for CAMI and how the code was up to 2017.
 * but, here's the twist, it turns out that this code doesn't remove any SNPs. 
 * kinda embarassing. metagenomes don't have SNPs so that's why i never noticed so far.
 * (maybe should have noticed that my human genomes assemblies were so bad)
 *
 * let me write the analysis of spades that i should have written earlier:
 *
 * in bulge_removed.hpp, 
 *
 * if bulgecov > 1000, fail (i didn't implement this)
 *
 * delta = CountMaxDifference(max_delta_, g_.length(e), max_relative_delta_); 
 * means that delta is set to max(bulgelen*0.1, 3)
 *
 * so we're looking for alternative paths which are of length [bulgelen-delta;bulgelen+delta]
 *
 */ 
template<typename GraphType, typename Node, typename Edge>
unsigned long Simplifications<GraphType,Node,Edge>::removeBulges()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int maxBulgeLength = std::max((unsigned int)((double)k * _bulgeLen_kMult), (unsigned int)(k + _bulgeLen_kAdd)); // SPAdes, exactly
    unsigned int backtrackingLimit = k+_bulgeAltPath_kAdd;//maxBulgeLength; // arbitrary, but if too high it will take much time; // with unitigs, no reason that it has to depend on k, but for some reason, setting it to just "k" doesnt remove nearly as many bulges as k=20. todo investigate that someday.
    unsigned int altPathCovMult = _bulgeAltPath_covMult;

    // stats
    //
    unsigned long nbBulgesCandidates = 0;
    unsigned long nbBulgesRemoved = 0;
    unsigned long nbSimplePaths = 0;
    unsigned long nbLongSimplePaths = 0;
    unsigned long nbShortSimplePaths = 0;
    unsigned long nbTopologicalBulges = 0;
    unsigned long nbFirstNodeDeleted = 0;
    unsigned long nbFirstNodeGraphDeleted = 0;
    unsigned long nbNoAltPathBulgesLoop = 0, nbNoAltPathBulgesDepth = 0,nbNoAltPathBulgesDeadend = 0;
    unsigned long nbBadCovBulges = 0;

    unsigned long timeAll = 0, timePathFinding = 0, timeFailedPathFinding = 0, timeLongestFailure = 0,
                  timeSimplePath = 0, timeDelete = 0, timePost = 0, timeVarious = 0;

    unsigned long longestFailureDepth = 0;

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, simplprogressFormat2, ++_nbBulgeRemovalPasses);
    ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem> *itNode; 
    if (_firstNodeIteration )
    {
        itNode = new ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem>(_graph.GraphType::iterator(), buffer, _verbose);
        if (_verbose)
            std::cout << "iterating on " << itNode->size() << " nodes" << std::endl;
    }
    else
    {
        itNode = new ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem>(_graph.GraphType::iteratorCachedNodes(), buffer, _verbose);
        if (_verbose)
            std::cout << "iterating on " << itNode->size() << " cached nodes" << std::endl;
    }


    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);

    NodesDeleter<Node,Edge,GraphType> nodesDeleter(_graph, nbNodes, _nbCores, _verbose);

#ifdef SIMPLIFICATION_LAMBDAS 
    dispatcher.iterate (itNode, [&] (Node& node) {
#else
    for (itNode->first(); !itNode->isDone(); itNode->next())
    {
        Node& node = itNode->item();
#endif
        TIME(auto start_thread_t=get_wtime());

        if (_graph.isNodeDeleted(node)) { 
#ifdef SIMPLIFICATION_LAMBDAS
            return; 
#else
            continue;
#endif
        }

      TIME(auto start_various_overhead_t=get_wtime());

          unsigned inDegree = _graph.indegree(node), outDegree = _graph.outdegree(node);

      TIME(auto end_various_overhead_t=get_wtime());
      TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));



      // need to search in both directions
      for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
      {
         if ((outDegree >= 2 && dir == DIR_OUTCOMING) || (inDegree >= 2 && dir == DIR_INCOMING))
         {
            __sync_fetch_and_add(&nbBulgesCandidates,1);

            TIME(auto start_various_overhead_t=get_wtime());
    
                /** We follow the outgoing simple paths to get their length and last neighbor */
                GraphVector<Edge> neighbors = _graph.neighborsEdge(node, dir);
                
                DEBUG_BULGES(DebugBR debugBR(_graph.toString(node), neighbors.size()););


            TIME(auto end_various_overhead_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

            // do everying for each possible short simple path that is neighbor of that node
            assert(neighbors.size() > 1);
            for (unsigned int i = 0; i < neighbors.size(); i++)
            {
                if (node == (neighbors[i].to)) // node being it's own neighbors, is rare, but let's avoid it
                    continue;
            
                TIME(auto start_various_overhead_t=get_wtime());

                    if (_graph.isNodeDeleted(neighbors[i].to)) { 
                         __sync_fetch_and_add(&nbFirstNodeGraphDeleted, 1);
                        continue;}

                TIME(auto end_various_overhead_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

                /* that's the putative bulge*/
                TIME(auto start_simplepath_t=get_wtime());
                Node&     simplePathStart = neighbors[i].to;
                Direction simplePathDir   = dir;
                unsigned int pathLen = _graph.simplePathLength(simplePathStart,simplePathDir);
                TIME(auto end_simplepath_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeSimplePath, diff_wtime(start_simplepath_t,end_simplepath_t)));
                __sync_fetch_and_add(&nbSimplePaths, 1);

                bool isShort = true;

                if (k + pathLen > maxBulgeLength) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                {
                    __sync_fetch_and_add(&nbLongSimplePaths, 1);
                    isShort = false;
                }

                if (!isShort || pathLen == 0) // can't do much if it's pathLen=0, we don't support edge removal, only node removal
                    continue;
                
                __sync_fetch_and_add(&nbShortSimplePaths, 1);

                TIME(start_various_overhead_t=get_wtime());

                
                    Node lastNode = _graph.simplePathLastNode(simplePathStart,simplePathDir);
                    GraphVector<Edge> outneighbors = _graph.neighborsEdge(lastNode, dir);
    
                    if (outneighbors.size() == 0) // might still be a tip, unremoved for some reason
                        continue;
   
                    // so here is a hidden assumption: maybe outneighbors is of size more than 1, we used to care about the first noed after.
                    // i could decide to enforce bulge popping only if the outneighbor has size 1. 
                    //if (outneighbors.size() != 1) 
                    //    continue;
                    //  but i'm decided to do without for now. TODO: explore all the end nodes, not just the first once

                    Node endNode = outneighbors[0].to;
    
                    // at this point, the last node in "nodes" is the last node of a potential Bulge path, and endNode is hopefully a branching node right after.
                    // check if it's connected to something that has in-branching. 
                    bool isDoublyConnected = (dir==DIR_OUTCOMING && _graph.indegree(endNode) > 1) || (dir==DIR_INCOMING && _graph.outdegree(endNode) > 1);
    
                    bool isTopologicalBulge = isDoublyConnected;
    
                    DEBUG_BULGES(debugBR.infoNeighbor(i, _graph.toString(neighbors[i].to), DIR2STR(dir), _graph.toString(lastNode), pathLen, outneighbors.size(), _graph.toString(endNode), isTopologicalBulge););

                TIME(end_various_overhead_t=get_wtime());
                TIME(__sync_fetch_and_add(&timeVarious, diff_wtime(start_various_overhead_t,end_various_overhead_t)));

                if (!isTopologicalBulge)
                    continue;
                    
                //cout << "endnode has indegree/outdegree: " <<_graph.indegree(endNode) << "/" << _graph.outdegree(endNode) <<  " and the node before has indegree/outdegree: " <<_graph.indegree(nodes.back()) << "/" << _graph.outdegree(nodes.back()) << endl;

                __sync_fetch_and_add(&nbTopologicalBulges, 1);

                unsigned int maxlen = std::max((unsigned int)(pathLen * 1.1),(unsigned int) (pathLen + 3)); // following SPAdes
                double mean_abundance_most_covered;
                int success;
                Node startNode = node;

                TIME(auto start_pathfinding_t=get_wtime());

                Path_t<Node>  heuristic_p_most; // actually won't be used.. (it's just for debug) so would be nice to get rid of it someday, but i don't want to deal with pointers.

                /* startNode is branching, because we want to find alternative paths, except the one that go through (neighbors[i].to)*/
                this->heuristic_most_covered_path(dir, startNode, endNode, 
                            maxlen, 
                            success, mean_abundance_most_covered,
                            heuristic_p_most,
                            backtrackingLimit, // avoid too much backtracking
                            &(neighbors[i].to), // avoid that node
                            true, // most covered path
                            false // kmer version (true), unitigs version (false) 
                            );

                TIME(auto end_pathfinding_t=get_wtime());
                TIME(__sync_fetch_and_add(&timePathFinding, diff_wtime(start_pathfinding_t,end_pathfinding_t)));

                if (success != 1)
                {
                    DEBUG_BULGES(cout << "HMCP failed: " << hmcpstatus2ascii(success) << endl);
                    TIME(__sync_fetch_and_add(&timeFailedPathFinding, diff_wtime(start_pathfinding_t,end_pathfinding_t)));
                    TIME(if (diff_wtime(start_pathfinding_t,end_pathfinding_t) > timeLongestFailure) { timeLongestFailure = diff_wtime(start_pathfinding_t,end_pathfinding_t); });
                    longestFailureDepth = maxlen;

                    if (success == HMCP_LOOP)
                        __sync_fetch_and_add(&nbNoAltPathBulgesLoop, 1);
                    if (success == HMCP_MAX_DEPTH)
                         __sync_fetch_and_add(&nbNoAltPathBulgesDepth, 1);
                    if (success == HMCP_DIDNT_FIND_END )
                         __sync_fetch_and_add(&nbNoAltPathBulgesDeadend, 1);
                    continue;
                }

                TIME(auto start_post_t=get_wtime());

                    bool debug_hmcp = false;
                    if (debug_hmcp)
                    {

                        double mean_abundance_least_covered;
                        Path_t<Node> heuristic_p_least, heuristic_p_most;
                        this->heuristic_most_covered_path(dir, startNode, endNode, maxlen, success, mean_abundance_most_covered,  heuristic_p_most,  backtrackingLimit, &(neighbors[i].to),  true, true /* old version */);
                        this->heuristic_most_covered_path(dir, startNode, endNode, maxlen, success, mean_abundance_least_covered, heuristic_p_least, backtrackingLimit, &(neighbors[i].to), false,  true /* old version */);
                        cout << "alternative path is:  "<< this->path2string(dir, heuristic_p_most, endNode)<< " abundance: "<< mean_abundance_most_covered <<endl;
                        DEBUG_BULGES(cout << endl << "alternative least is: "<< this->path2string(dir, heuristic_p_least, endNode)<< " abundance: "<< mean_abundance_least_covered <<endl);
                    }

                    double simplePathCoverage = _graph.simplePathMeanAbundance(simplePathStart, simplePathDir);
    
                    bool isBulge =  simplePathCoverage <=  mean_abundance_most_covered * altPathCovMult /*typically 1.1 in genome assembly, SPAdes*/;
    
                    DEBUG_BULGES(cout << "bulge coverages: " << simplePathCoverage << " (path: " << _graph.toString(simplePathStart) << ") vs most covered:" <<  mean_abundance_most_covered  << endl);
    
                    if (!isBulge)
                    {
                        __sync_fetch_and_add(&nbBadCovBulges, 1);
                        DEBUG_BULGES(cout << "not a bulge due to coverage criterion" << endl);

                        TIME(auto end_post_t=get_wtime());
                        TIME(__sync_fetch_and_add(&timePost, diff_wtime(start_post_t,end_post_t)));
                        continue;
                    }

                    if (nodesDeleter.get(simplePathStart))
                    {
                        // not double-counting that delete
                    }
                    else
                    {
                        // delete the bulge
                        //
                        DEBUG_BULGES(cout << endl << "BULGE of length " << pathLen << " FOUND: " <<  _graph.toString (simplePathStart) << endl);
                        _graph.simplePathDelete(simplePathStart, simplePathDir, nodesDeleter);

                        __sync_fetch_and_add(&nbBulgesRemoved, 1);
                    }

                TIME(auto end_post_t=get_wtime());
                TIME(__sync_fetch_and_add(&timePost, diff_wtime(start_post_t,end_post_t)));

                break; // quite important to break here: don't try to remove the other neighbor (which might also satisfy the bulge condition)

            } // for neighbors

            DEBUG_BULGES(debugBR.draw(););
        } // if outdegree
      } // for direction
        TIME(auto end_thread_t=get_wtime());
        TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t)));
    
#ifdef SIMPLIFICATION_LAMBDAS
            }); // parallel
#else
}
#endif


    // now delete all nodes, in parallel
    TIME(auto start_nodedelete_t=get_wtime());

    nodesDeleter.flush();

    TIME(auto end_nodedelete_t=get_wtime());
    TIME(__sync_fetch_and_add(&timeDelete, diff_wtime(start_nodedelete_t,end_nodedelete_t)));

    if (_verbose)
    {
        cout << nbBulgesRemoved << " bulges removed. " << endl <<
            nbSimplePaths << "/" << nbLongSimplePaths << "+" <<nbShortSimplePaths << " any=long+short simple path examined across all threads, among them " <<
            nbTopologicalBulges << " topological bulges, " << nbFirstNodeDeleted << "+" << nbFirstNodeGraphDeleted << " were first-node duplicates." << endl;
        cout << nbBulgesCandidates << " bulges candidates passed degree check. " << endl;
        cout << nbNoAltPathBulgesDepth << "+" << nbNoAltPathBulgesLoop << "+" << nbNoAltPathBulgesDeadend << " without alt. path (complex+loop+noend), " 
            << nbBadCovBulges << " didn't satisfy cov. criterion." << endl;
    }

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);

    if (_verbose)
    {
        TIME(cout << "Bulges timings: " << timeAll / unit << " CPUsecs total."<< endl);
        TIME(cout << "                " << timeSimplePath / unit << " CPUsecs simple path traversal." << endl);
        TIME(cout << "                " << timePathFinding / unit << "(/" << timePathFinding / unit << ") CPUsecs path-finding(/failed). Longest: " << timeLongestFailure / (unit/1000) << " CPUmillisecs (depth " << longestFailureDepth << ")." << endl);
        TIME(cout << "                " << timePost / unit << " CPUsecs topological bulge processing, " << endl);
        TIME(cout << "                " << timeDelete / unit << " CPUsecs nodes deletion." << endl);
        TIME(cout << "                " << timeVarious / unit << " CPUsecs various overhead." << endl);
    }

    return nbBulgesRemoved;
}



/* this is an algorithm that simplifies the graph, again taken from SPAdes
 * it removes so-called "erroneous connections" which are Z-shaped motifs in the graph
 *
 *
 * Again, let's see what spades 3.5 does.
 *erroneous connection remover is:

 RemoveLowCoverageEdges
 calls config then calls:
  omnigraph::EdgeRemovingAlgorithm which is same as tc

  to_ec_lb in ../src/debruijn/simplification/simplification_settings.hpp
  is exactly like tip clipping but with a different length (2 * (max tip length with coeff 5) - 1, so that's exactly 10*min(k,readlen/2) - 1)

  icb is something that removes edges with coverage = (cov_bound * iter / nb_iters )
  where cov_bound is same as in cb

  and it's asking for AlternativesPresenceCondition:

    o   o                           o-->-O
   /     \                             /
   O-->---O    drawn differently:     /
                                    O-->-o

  so anyway, since we're not computing the coverage model like SPAdes does, I'm going to use RCTC 4 (found that RCTC 2 gives too many misassemblies)

*/
template<typename GraphType, typename Node, typename Edge>
unsigned long Simplifications<GraphType,Node,Edge>::removeErroneousConnections()
{
    unsigned int k = _graph.getKmerSize();
    unsigned int maxECLength = (unsigned int)((float)k * _ecLen_kMult) ;  // SPAdes mode 
    double RCTCcutoff = _ecRCTCcutoff;

    unsigned long nbSimplePaths = 0;
    unsigned long nbLongSimplePaths = 0;
    unsigned long nbShortSimplePaths = 0;
    unsigned long nbTopologicalEC = 0;
    unsigned long nbECRemoved = 0;
    unsigned long nbECCandidates = 0;
    unsigned long timeDelete = 0, timeProcessing = 0, timeCoverage = 0;
    unsigned long timeAll = 0, timeSimplePath = 0;

    /** We get an iterator over all nodes . */
    char buffer[128];
    sprintf(buffer, simplprogressFormat3, ++_nbECRemovalPasses);
    ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem> *itNode; 
    if (_firstNodeIteration )
    {
        itNode = new ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem>(_graph.GraphType::iterator(), buffer, _verbose);
        if (_verbose)
            std::cout << "iterating on " << itNode->size() << " nodes on disk" << std::endl;
    }
    else
    {
        itNode = new ProgressGraphIteratorTemplate<Node,ProgressTimerAndSystem>(_graph.GraphType::iteratorCachedNodes(), buffer, _verbose);
        if (_verbose)
            std::cout << "iterating on " << itNode->size() << " cached nodes" << std::endl;
    }

    // parallel stuff: create a dispatcher ; support atomic operations
    Dispatcher dispatcher (_nbCores);

    // parallel stuff
    NodesDeleter<Node,Edge,GraphType> nodesDeleter(_graph, nbNodes, _nbCores, _verbose);

#ifdef SIMPLIFICATION_LAMBDAS 
    dispatcher.iterate (itNode, [&] (Node& node) {
#else
    for (itNode->first(); !itNode->isDone(); itNode->next())
    {
        Node& node = itNode->item();
#endif
            TIME(auto start_thread_t=get_wtime());

            if (_graph.isNodeDeleted(node)) { 
#ifdef SIMPLIFICATION_LAMBDAS
            return; 
#else
            continue;
#endif
            }

            unsigned inDegree = _graph.indegree(node), outDegree = _graph.outdegree(node);

            /* ec nodes have out/in degree of 1 or more on one side, and of 2 or more on the other */
            if (!((inDegree >= 1 && outDegree > 1 ) || (inDegree > 1 && outDegree >=1 )))
            {
#ifdef SIMPLIFICATION_LAMBDAS
            return; 
#else
            continue;
#endif
            }

                /* at this point, "node" is a node such that:
                 *
                 *            /                  \
                 *           /                    \
                 *---- [node] -----     or  -------[node]-----
                 * so it's the potential node right before an EC. 
                 */

            // need to search in both directions
            for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
            {
                if ((outDegree >= 2 && dir == DIR_OUTCOMING) || (inDegree >= 2 && dir == DIR_INCOMING))
                {  
                    DEBUG_EC(cout << endl << "putative EC node: " << _graph.toString (node) << endl);
                    __sync_fetch_and_add(&nbECCandidates,1);

                    /** We follow the outcoming simple paths 
                     * (so, if it's outdegree 2, we follow them to get their length and last neighbor */
                    GraphVector<Edge> neighbors = _graph.neighborsEdge(node, dir);

                    assert(neighbors.size() > 1);

                    // do everying for each possible short simple path that is neighbor of that node
                    for (unsigned int i = 0; i < neighbors.size(); i++)
                    {
                        if (_graph.isNodeDeleted(neighbors[i].to)) { 
                            continue;}

                        /* at this point, neighbors[i].to is the first node of an EC, maybe. */
                         
                        /* explore the simple path from that node */
                        TIME(auto start_simplepath_t=get_wtime());
                        Node&     simplePathStart = neighbors[i].to;
                        Direction simplePathDir   = dir;
                        unsigned int pathLen = _graph.simplePathLength(simplePathStart,simplePathDir);
                        TIME(auto end_simplepath_t=get_wtime());
                        TIME(__sync_fetch_and_add(&timeSimplePath, diff_wtime(start_simplepath_t,end_simplepath_t)));

                        DEBUG_EC(cout << endl << "neighbors " << i+1 << "/" << neighbors.size() << " from: " << _graph.toString (neighbors[i].to) << " dir: " << DIR2STR(dir) << endl);
                        bool isShort = true;
                        
                        if (k + pathLen > maxECLength) // "k +" is to take into account that's we're actually traversing a path of extensions from "node"
                        {
                            __sync_fetch_and_add(&nbLongSimplePaths, 1);
                            isShort = false;
                        }

                        __sync_fetch_and_add(&nbSimplePaths, 1);

                        if ((!isShort) || pathLen == 0) // can't do much if it's pathLen=0, we don't support edge removal, only node removal
                        {
                            DEBUG_EC(cout << "direction: " << DIR2STR(dir) << ", not an EC: isShort: " << isShort << " pathLen: " << pathLen << endl);
                            continue;
                        }

                        __sync_fetch_and_add(&nbShortSimplePaths, 1);
            
                        Node lastNode           = _graph.simplePathLastNode     (simplePathStart,simplePathDir);
                        double pathMeanAbundance = _graph.simplePathMeanAbundance(simplePathStart,simplePathDir);

                        GraphVector<Edge> outneighbors = _graph.neighborsEdge(lastNode, dir);
                        DEBUG_EC(cout << "last simple path node: "<< _graph.toString(lastNode) << " has " << outneighbors.size() << " outneighbors; mean abundance: " << pathMeanAbundance << endl);

                        if (outneighbors.size() == 0) // might still be a tip, unremoved for some reason
                            continue;

                        Node endNode = outneighbors[0].to;
                        DEBUG_EC(cout << "endNode: " << _graph.toString(endNode) << endl);

                        // at this point, the last node in "nodes" is the last node of a potential EC, and endNode is hopefully a branching node right after.
                        // check if it's connected to something that has in-branching and also an out neighbor. 
                        bool isDoublyConnected = (dir==DIR_OUTCOMING && _graph.indegree(endNode) > 1 && _graph.outdegree(endNode) >= 1) || \
                                                 (dir==DIR_INCOMING && _graph.outdegree(endNode) > 1 && _graph.indegree(endNode) >= 1 );


                        bool isTopologicalEC = isDoublyConnected;

                        DEBUG_EC(cout << "direction: " << DIR2STR(dir) << ", pathlen: " << pathLen << " last node neighbors size: " << _graph.neighborsEdge(lastNode).size() << " indegree outdegree: " <<_graph.indegree(node) << " " << _graph.outdegree(node) << " isDoublyConnected: " << isDoublyConnected << " isTopoEC: " << isTopologicalEC << endl);


                        if (isTopologicalEC)
                        {
                            __sync_fetch_and_add(&nbTopologicalEC,1);
                            TIME(auto start_ec_coverage_t=get_wtime());
                            bool isRCTC = this->satisfyRCTC(pathMeanAbundance, lastNode, RCTCcutoff, dir);

                            isRCTC |= this->satisfyRCTC(pathMeanAbundance, neighbors[i].to, RCTCcutoff, reverse(dir)); // also check in the other direction
                            // TODO think hard, is it a |= or a &=? FIXME for potential misassemblies
                            
                            TIME(auto end_ec_coverage_t=get_wtime());
                            TIME(__sync_fetch_and_add(&timeCoverage, diff_wtime(start_ec_coverage_t,end_ec_coverage_t)));

                            bool isEC = isRCTC;
                        
                            DEBUG_EC(cout << "isRCTC:" << isRCTC << endl);

                            TIME(auto start_ec_processing_t=get_wtime());

                            if (isEC)
                            {
                                if (nodesDeleter.get(simplePathStart))
                                {
                                    // not double-counting that delete
                                }
                                else
                                {
                                    // delete it
                                    //
                                    _graph.simplePathDelete(simplePathStart, simplePathDir, nodesDeleter);
                                    DEBUG_EC(cout << endl << "EC of length " << pathLen << " FOUND: " <<  _graph.toString (node) << endl);

                                    __sync_fetch_and_add(&nbECRemoved, 1);
                                }

                            }
                            TIME(auto end_ec_processing_t=get_wtime());
                            TIME(__sync_fetch_and_add(&timeProcessing, diff_wtime(start_ec_processing_t,end_ec_processing_t)));
                        }
                    }
                }
            }
            TIME(auto end_thread_t=get_wtime());
            TIME(__sync_fetch_and_add(&timeAll, diff_wtime(start_thread_t,end_thread_t)));

#ifdef SIMPLIFICATION_LAMBDAS
            }); // parallel
#else
}
#endif

    TIME(auto start_nodesdel_t=get_wtime());

    nodesDeleter.flush();

    TIME(auto end_nodesdel_t=get_wtime()); 
    TIME(__sync_fetch_and_add(&timeDelete, diff_wtime(start_nodesdel_t,end_nodesdel_t)));

    if (_verbose)
    {
        cout << nbECRemoved << " erroneous connections removed. " << endl;
        cout << nbECCandidates << " EC candidates passed degree check. " << endl;
        cout << nbSimplePaths << "/" << nbLongSimplePaths << "+" <<nbShortSimplePaths << " any=long+short simple path examined across all threads" << endl;
        cout << nbTopologicalEC << " topological ECs. " << endl;
        double unit = 1000000000;
        cout.setf(ios_base::fixed);
        cout.precision(1);
        TIME(cout << "EC Timings: " << timeAll / unit << " CPUsecs total, including"<< endl);
        TIME(cout << "                " << timeSimplePath / unit << " CPUsecs EC simple paths" << endl);
        TIME(cout << "                " << timeCoverage / unit << " CPUsecs EC coverage test" << endl);
        TIME(cout << "                " << timeProcessing / unit << " CPUsecs EC processing" << endl);
        TIME(cout << "Nodes deletion: " << timeDelete / unit << " CPUsecs." << endl);
    }

    return nbECRemoved;
}

// instantiation
template class Simplifications<Graph, Node, Edge>; 

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
