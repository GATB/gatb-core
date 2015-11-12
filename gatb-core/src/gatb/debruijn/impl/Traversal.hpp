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

#ifndef _GATB_TOOLS_TRAVERSAL_HPP_
#define _GATB_TOOLS_TRAVERSAL_HPP_

#include <gatb/debruijn/impl/Terminator.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <set>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

// some stats
struct TraversalStats
{
    // MonumentTraversal
    long ended_traversals;
    long couldnt_find_all_consensuses;
    long couldnt_validate_consensuses;
    long couldnt_traverse_bubble_breadth;
    long couldnt_traverse_bubble_depth;
    long couldnt_because_marked_kmer;
    long couldnt_inbranching_depth;
    long couldnt_inbranching_breadth;
    long couldnt_inbranching_other;
    long couldnt_find_extension;
    //find_all_consensuses failures:
    long couldnt_consensus_negative_depth;
    long couldnt_consensus_amount;
    long couldnt_consensus_loop;
    //validate_consensus potentials failures:
    long couldnt_validate_bubble_mean_depth;
    long couldnt_validate_bubble_stdev;
    long couldnt_validate_bubble_deadend;
    long couldnt_validate_bubble_identity;
    long couldnt_validate_bubble_long_chosen;


    //SimplepathTraversal
    long couldnt_no_extension;
    long couldnt_inbranching;
    long couldnt_outbranching;
};

/********************************************************************************/

/** \brief Class that traverse nodes of a GraphTemplate<Node,Edge,GraphDataVariant>
 *
 * The Traversal class looks for path in a graph according to several criteria. This
 * is done through the 'traverse' method. As a result, one gets a Path object that
 * holds the traversed path information.
 *
 * There are two kinds of traversal:
 *  - \ref tools::misc::TRAVERSAL_UNITIG : one gets only simple path in the graph
 *  - \ref tools::misc::TRAVERSAL_CONTIG : one gets more complex path
 *
 * A factory method \ref create is available and should be used by end users to instantiate
 * this class.
 *
 * This class is abstract since it doesn't implement avance. It is subclassed according
 * to the wanted kind of traversal.
 *
 * Example of use : we create a fake bank, build its graph and traverse the graph:
 * \snippet traversal2.cpp  snippet1_traversal
 *
 */
template <typename Node, typename Edge, typename GraphDataVariant>
class TraversalTemplate: public system::SmartPointer
{
public:

    /** Factory method that creates an instance of Traversal
     * \param[in] type : kind of traversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] max_len : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    static TraversalTemplate<Node,Edge,GraphDataVariant>* create (
        tools::misc::TraversalKind  type,
        const GraphTemplate<Node,Edge,GraphDataVariant>&                graph,
        TerminatorTemplate<Node,Edge,GraphDataVariant>&                         terminator,
        int                         max_len     = defaultMaxLen,
        int                         max_depth   = defaultMaxDepth,
        int                         max_breadth = defaultMaxBreadth
    );

    /** Factory method that creates an instance of Traversal
     * \param[in] type : type of traversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] max_len : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    static TraversalTemplate<Node,Edge,GraphDataVariant>* create (
        const std::string&  type,
        const GraphTemplate<Node,Edge,GraphDataVariant>&        graph,
        TerminatorTemplate<Node,Edge,GraphDataVariant>&                 terminator,
        int                 max_len     = defaultMaxLen,
        int                 max_depth   = defaultMaxDepth,
        int                 max_breadth = defaultMaxBreadth
    );

    /** Get the name of the traversal.
     * \return traversal name.  */
    virtual std::string getName() const = 0;

    /** Traverse the graph starting from one node for a given direction. As a result,
     * one gets a Path object
     * \param[in] node : starting node of the traversal.
     * \param[in] dir :  direction of the traversal
     * \param[out] resulting_sequence : path of the traversal */
    int traverse (Node& node, Direction dir, Path_t<Node>& resulting_sequence)
    {
        Node endingNode;  return traverse (node, endingNode, dir, resulting_sequence);
    }

    /** Traverse the graph starting from one node for a given direction. As a result,
     * one gets a Path object. We get also the last traversed node.
     * \param[in] startingNode : starting node of the traversal.
     * \param[out] endingNode : starting node of the traversal.
     * \param[in] dir :  direction of the traversal
     * \param[out] resulting_sequence : path of the traversal */
    int traverse (Node& startingNode, Node& endingNode, Direction dir, Path_t<Node>& resulting_sequence);

    /** Get the maximum allowed depth
     * \return maximum depth */
    unsigned int getMaxDepth() const  { return max_depth; };

    /** Get the maximum allowed breadth
     * \return maximum breadth */
    unsigned int getMaxBreadth () const  { return max_breadth; };

    static const int defaultMaxLen     = 10*1000*1000;
    static const int defaultMaxDepth   = 500;
    static const int defaultMaxBreadth = 20;

    TraversalStats final_stats, stats;
    void commit_stats() { final_stats = stats; }; // save current stats into final stats
    void revert_stats() { stats = final_stats; }; // discard changes in stats (because contig was discarded)

    /** Compute a global alignment between two path. NOTE: could be moved to Path class.
     * \param[in] a : first path
     * \param[in] b : second path. */
    static float needleman_wunch (const Path_t<Node>& a, const Path_t<Node>& b);

    /** Get the bubbles found during traversal. One bubble is defined by the [begin,end] positions in the
     * path.
     * \return vector of positions ranges. */
    const std::vector <std::pair<int, int> >& getBubbles()  const { return bubbles_positions; }

    bool deadend;

protected:

    /** */
    TraversalTemplate (
        const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
        TerminatorTemplate<Node,Edge,GraphDataVariant>         & terminator,
        int maxlen,
        int max_depth,
        int max_breadth
    );

    const GraphTemplate<Node,Edge,GraphDataVariant>& graph;
    TerminatorTemplate<Node,Edge,GraphDataVariant>&          terminator;

    int maxlen;
    int max_depth;
    int max_breadth;

    virtual char avance (Node& node, Direction dir, bool first_extension, Path_t<Node>& path, Node& previousNode) = 0;

    void mark_extensions (std::set<Node>& extensions_to_mark);

    // record the start/end positions of traversed bubbles (only from the latest traverse() call)
    std::vector <std::pair<int, int> > bubbles_positions;
};

/********************************************************************************/

/** \brief Null implementation of Traversal.
 *
 * This class returns empty Path as a result of traverse.
 */
template <typename Node, typename Edge, typename GraphDataVariant>
class NullTraversalTemplate: public TraversalTemplate<Node,Edge,GraphDataVariant>
{
public:

    /** Factory method that creates an instance of NullTraversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] maxlen : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    NullTraversalTemplate (
        const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
        TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator,
        int maxlen      = NullTraversalTemplate::defaultMaxLen,
        int max_depth   = NullTraversalTemplate::defaultMaxDepth,
        int max_breadth = NullTraversalTemplate::defaultMaxBreadth
    ) : TraversalTemplate<Node,Edge,GraphDataVariant> (graph, terminator, maxlen, max_depth, max_breadth) {}

    /** Get the name of the traversal
     * \return the name */
    std::string getName() const  { return tools::misc::toString(tools::misc::TRAVERSAL_NONE); }

private:

    char avance (Node& node, Direction dir, bool first_extension, Path_t<Node>& path, Node& previousNode) { return 0; }
};

/********************************************************************************/

/** \brief Implementation of Traversal that produces unitigs.
 */
template <typename Node, typename Edge, typename GraphDataVariant>
class SimplePathsTraversalTemplate: public TraversalTemplate<Node,Edge,GraphDataVariant>
{
public:

    /** Factory method that creates an instance of SimplePathsTraversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] maxlen : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    SimplePathsTraversalTemplate (
        const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
        TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator,
        int maxlen      = SimplePathsTraversalTemplate::defaultMaxLen,
        int max_depth   = SimplePathsTraversalTemplate::defaultMaxDepth,
        int max_breadth = SimplePathsTraversalTemplate::defaultMaxBreadth
    );

    /** Get the name of the traversal
     * \return the name */
    std::string getName() const  { return tools::misc::toString(tools::misc::TRAVERSAL_UNITIG); }

private:

    /* Implementation of the virtual method. */
    char avance (Node& node, Direction dir, bool first_extension, Path_t<Node>& path, Node& previousNode);
};

/********************************************************************************/

/** \brief Implementation of Traversal that produces contigs.
 */
template <typename Node, typename Edge, typename GraphDataVariant>
class MonumentTraversalTemplate: public TraversalTemplate<Node,Edge,GraphDataVariant>
{
public:

    /** Factory method that creates an instance of MonumentTraversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] maxlen : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    MonumentTraversalTemplate (
        const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
        TerminatorTemplate<Node,Edge,GraphDataVariant>&          terminator,
        int maxlen      = MonumentTraversalTemplate::defaultMaxLen,
        int max_depth   = MonumentTraversalTemplate::defaultMaxDepth,
        int max_breadth = MonumentTraversalTemplate::defaultMaxBreadth
    );

    /** Get the name of the traversal
     * \return the name */
    std::string getName() const  { return tools::misc::toString(tools::misc::TRAVERSAL_CONTIG); }

    bool explore_branching (
        Node& node,
        Direction dir,
        Path_t<Node>& consensus,
        Node& previousNode,
        std::set<Node>& all_involved_extensions
    );


    // those two used to be private, but I need them in graph Simplifications for now (until explore_branching gets templated or any way we can choose Frontline type)
    bool validate_consensuses (std::set<Path_t<Node> >& consensuses, Path_t<Node>& consensus);

        std::set<Path_t<Node> > all_consensuses_between (
        Direction    dir,
        Node& startNode,
        Node& endNode,
        int traversal_depth,
        bool &success
    );


private:

    /* Implementation of the virtual method. */
    char avance (Node& node, Direction dir, bool first_extension, Path_t<Node>& path, Node& previousNode);

    bool explore_branching (
        Node& node,
        Direction dir,
        Path_t<Node>& consensus,
        Node& previousNode
    );

    int find_end_of_branching (
        Direction dir,
        Node& startingNode,
        Node& endNode,
        Node& previousNode,
        std::set<Node>& all_involved_extensions
    );
 
    std::set<Path_t<Node> > all_consensuses_between (
        Direction    dir,
        Node& startNode,
        Node& endNode,
        int traversal_depth,
        std::set<typename Node::Value> usedNode,
        Path_t<Node> current_consensus,
        bool& success
    );
   
    bool all_consensuses_almost_identical (std::set<Path_t<Node> >& consensuses);

    void mark_extensions (std::set<Node>& extensions_to_mark);

    Path_t<Node> most_abundant_consensus(std::set<Path_t<Node> >& consensuses);

    static const int consensuses_identity = 80; // traversing bubble if paths are all pair-wise identical by 80% 
    //(used to be > 90% in legacy minia) // by legacy minia i mean minia 1 and minia 2 up to the assembly algo rewrite in may 2015
};

/* typedef for compatibility with all existing GATB tools */

typedef TraversalTemplate<Node, Edge, GraphDataVariant> Traversal; 
typedef MonumentTraversalTemplate<Node, Edge, GraphDataVariant> MonumentTraversal; 
typedef NullTraversalTemplate<Node, Edge, GraphDataVariant> NullTraversal; 
typedef SimplePathsTraversalTemplate<Node, Edge, GraphDataVariant> SimplePathsTraversal;



/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_TRAVERSAL_HPP_ */

