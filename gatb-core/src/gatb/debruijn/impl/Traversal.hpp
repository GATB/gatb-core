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
};

/********************************************************************************/

/** \brief Class that traverse nodes of a Graph
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
class Traversal : public system::SmartPointer
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
    static Traversal* create (
        tools::misc::TraversalKind  type,
        const Graph&                graph,
        Terminator&                 terminator,
        int                         max_len     = Traversal::defaultMaxLen,
        int                         max_depth   = Traversal::defaultMaxDepth,
        int                         max_breadth = Traversal::defaultMaxBreadth
    );

    /** Factory method that creates an instance of Traversal
     * \param[in] type : type of traversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] max_len : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    static Traversal* create (
        const std::string&  type,
        const Graph&        graph,
        Terminator&         terminator,
        int                 max_len     = Traversal::defaultMaxLen,
        int                 max_depth   = Traversal::defaultMaxDepth,
        int                 max_breadth = Traversal::defaultMaxBreadth
    );

    /** Get the name of the traversal.
     * \return traversal name.  */
    virtual std::string getName() const = 0;

    /** Traverse the graph starting from one node for a given direction. As a result,
     * one gets a Path object
     * \param[in] node : starting node of the traversal.
     * \param[in] dir :  direction of the traversal
     * \param[out] resulting_sequence : path of the traversal */
    int traverse (const Node& node, Direction dir, Path& resulting_sequence)
    {
        Node endingNode;  return traverse (node, endingNode, dir, resulting_sequence);
    }

    /** Traverse the graph starting from one node for a given direction. As a result,
     * one gets a Path object. We get also the last traversed node.
     * \param[in] startingNode : starting node of the traversal.
     * \param[out] endingNode : starting node of the traversal.
     * \param[in] dir :  direction of the traversal
     * \param[out] resulting_sequence : path of the traversal */
    int traverse (const Node& startingNode, Node& endingNode, Direction dir, Path& resulting_sequence);

    /** Get the maximum allowed depth
     * \return maximum depth */
    int getMaxDepth() const  { return max_depth; };

    /** Get the maximum allowed breadth
     * \return maximum breadth */
    int getMaxBreadth () const  { return max_breadth; };

    static const int defaultMaxLen     = 10*1000*1000;
    static const int defaultMaxDepth   = 500;
    static const int defaultMaxBreadth = 20;

    TraversalStats final_stats, stats;
    void commit_stats() { final_stats = stats; }; // save current stats into final stats
    void revert_stats() { stats = final_stats; }; // discard changes in stats (because contig was discarded)

    /** Compute a global alignment between two path. NOTE: could be moved to Path class.
     * \param[in] a : first path
     * \param[in] b : second path. */
    static float needleman_wunch (const Path& a, const Path& b);

    /** Get the bubbles found during traversal. One bubble is defined by the [begin,end] positions in the
     * path.
     * \return vector of positions ranges. */
    const std::vector <std::pair<int, int> >& getBubbles()  const { return bubbles_positions; }

protected:

    /** */
    Traversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen,
        int max_depth,
        int max_breadth
    );

    const Graph& graph;
    Terminator&  terminator;

    int maxlen;
    int max_depth;
    int max_breadth;

    virtual char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode) = 0;

    void mark_extensions (std::set<Node>& extensions_to_mark);

    // record the start/end positions of traversed bubbles (only from the latest traverse() call)
    std::vector <std::pair<int, int> > bubbles_positions;
};

/********************************************************************************/

/** \brief Null implementation of Traversal.
 *
 * This class returns empty Path as a result of traverse.
 */
class NullTraversal: public Traversal
{
public:

    /** Factory method that creates an instance of NullTraversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] maxlen : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    NullTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = Traversal::defaultMaxLen,
        int max_depth   = Traversal::defaultMaxDepth,
        int max_breadth = Traversal::defaultMaxBreadth
    ) : Traversal (graph, terminator, maxlen, max_depth, max_breadth) {}

    /** Get the name of the traversal
     * \return the name */
    std::string getName() const  { return tools::misc::toString(tools::misc::TRAVERSAL_NONE); }

private:

    char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode) { return 0; }
};

/********************************************************************************/

/** \brief Implementation of Traversal that produces unitigs.
 */
class SimplePathsTraversal: public Traversal
{
public:

    /** Factory method that creates an instance of SimplePathsTraversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] maxlen : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    SimplePathsTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = Traversal::defaultMaxLen,
        int max_depth   = Traversal::defaultMaxDepth,
        int max_breadth = Traversal::defaultMaxBreadth
    );

    /** Get the name of the traversal
     * \return the name */
    std::string getName() const  { return tools::misc::toString(tools::misc::TRAVERSAL_UNITIG); }

private:

    /* Implementation of the virtual method. */
    char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode);
};

/********************************************************************************/

/** \brief Implementation of Traversal that produces contigs.
 */
class MonumentTraversal: public Traversal
{
public:

    /** Factory method that creates an instance of MonumentTraversal
     * \param[in] graph : graph object to be traversed
     * \param[in] terminator : object used to tag traversed nodes
     * \param[in] maxlen : maximum length of the traversal
     * \param[in] max_depth : maximum depth of the traversal
     * \param[in] max_breadth : maximum depth of the traversal
     */
    MonumentTraversal (
        const Graph& graph,
        Terminator& terminator,
        int maxlen      = Traversal::defaultMaxLen,
        int max_depth   = Traversal::defaultMaxDepth,
        int max_breadth = Traversal::defaultMaxBreadth
    );

    /** Get the name of the traversal
     * \return the name */
    std::string getName() const  { return tools::misc::toString(tools::misc::TRAVERSAL_CONTIG); }

    bool explore_branching (
        const Node& node,
        Direction dir,
        Path& consensus,
        const Node& previousNode,
        std::set<Node>& all_involved_extensions
    );

private:

    /* Implementation of the virtual method. */
    char avance (const Node& node, Direction dir, bool first_extension, Path& path, const Node& previousNode);

    bool explore_branching (
        const Node& node,
        Direction dir,
        Path& consensus,
        const Node& previousNode
    );


    int find_end_of_branching (
        Direction dir,
        const Node& startingNode,
        Node& endNode,
        const Node& previousNode,
        std::set<Node>& all_involved_extensions
    );

    std::set<Path> all_consensuses_between (
        Direction    dir,
        const Node& startNode,
        const Node& endNode,
        int traversal_depth,
        std::set<Node::Value> usedNode,
        Path current_consensus,
        bool& success
    );

    std::set<Path> all_consensuses_between (
        Direction    dir,
        const Node& startNode,
        const Node& endNode,
        int traversal_depth,
        bool &success
    );

    bool validate_consensuses (std::set<Path>& consensuses, Path& consensus);

    bool all_consensuses_almost_identical (std::set<Path>& consensuses);

    void mark_extensions (std::set<Node>& extensions_to_mark);

    Path most_abundant_consensus(std::set<Path>& consensuses);

    static const int consensuses_identity = 90; // traversing bubble if paths are all pair-wise identical by > 90%
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_TRAVERSAL_HPP_ */

