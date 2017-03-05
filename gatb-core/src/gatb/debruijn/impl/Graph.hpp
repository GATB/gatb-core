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

/** \file Graph.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Graph class
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HPP_

/********************************************************************************/
#include <vector>
#include <set>

#include <unordered_map>

#include <gatb/system/api/IThread.hpp> // for ISynchronizer
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/math/Integer.hpp>

#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/kmer/impl/MPHFAlgorithm.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <gatb/debruijn/impl/NodesDeleter.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for De Bruijn graph management. */
namespace debruijn  {
/** \brief Implementation package for De Bruijn graph management. */
namespace impl      {

/********************************************************************************/

enum Direction
{
    DIR_OUTCOMING = 1,
    DIR_INCOMING = 2,
    DIR_END = 3
};

inline std::string toString (Direction d)
{
    if (d==DIR_OUTCOMING)      { return std::string("OUT"); }
    else if (d==DIR_INCOMING)  { return std::string("IN");  }
    else { return std::string ("???"); }
}

inline Direction reverse (Direction dir)  { return dir==DIR_OUTCOMING ? DIR_INCOMING : DIR_OUTCOMING; }

// Quite ugly... should be improved...
#define foreach_direction(d)  for (Direction d=DIR_OUTCOMING; d<DIR_END; d = (Direction)((int)d + 1) )


/********************************************************************************
                        #     #  #######  ######   #######
                        ##    #  #     #  #     #  #
                        # #   #  #     #  #     #  #
                        #  #  #  #     #  #     #  #####
                        #   # #  #     #  #     #  #
                        #    ##  #     #  #     #  #
                        #     #  #######  ######   #######
********************************************************************************/

/** \brief Node structure representing a node in the De Bruijn graph
 *
 * The Node structure needs at least two attributes for uniquely defining a node in the DBG:
 *  - a kmer value (representing a node and its reverse complement)
 *  - a strand value (telling on which strand the kmer has to be interpreted)
 *
 *  A specific Node::Value type is defined for the kmer type. It may be used sometimes when one
 *  wants to perform some action on a node without the need to know the strand.
 *
 *  Although it is possible to instantiate "from scratch" a Node object, it is likely that Node
 *  objects are retrieved through the Graph class. Indeed, a Node object is likely to be useful
 *  if it belongs to a graph, so it makes sense that the Graph class is the main actor for providing
 *  Node objects.
 *
 *  The Node structure may be inherited by other structure, mainly for refining some characteristics
 *  for some kinds of nodes. For instance, a BranchingNode structure may be defined for nodes that have
 *  more than one incoming and/or outcoming neighbors. Such inherited structures could be understood as
 *  regular Node instance with a specific invariant. They can also be the material for specialization of
 *  some template methods of the Graph class (see for instance Graph::iterator template specialization
 *  for the BranchingNode structure).
 *
 *  The Node structure has also a 'abundance' attribute; it gives the occurrences number of the kmer in
 *  the initial set of reads. Note that this attribute has a correct value only when nodes are iterated 
 *  (from the solid kmers on disk). When Node objects are created via the Graph API (e.g. using
 *  as a result of buildNode, successors(), predecessors(), etc..), developers should not assume
 *  that the 'abundance' attribute is correctly set. Use the queryAbundance() function instead to get
 *  the abundance of a node.
 */

template <typename Value_t=tools::math::Integer>
struct Node_t
{
    /** Type for a kmer value. */
    typedef Value_t Value;

    /** Default constructor. */
    Node_t() : strand(kmer::STRAND_FORWARD), abundance(0), mphfIndex(0) , iterationRank(0) {}

    /** Constructor.
     * \param[in] kmer : kmer value. By default, it is the minimum value of the forward and revcomp value.
     * \param[in] strand : strand telling how to interpret the kmer value. Default value is forward by convention.
     * \param[in] abundance : abundance of the kmer. Default value is 0 if not set.
     */
    Node_t (const Node_t::Value& kmer, kmer::Strand strand=kmer::STRAND_FORWARD, u_int16_t abundance=0, u_int64_t mphfIndex = 0)
        : kmer(kmer), strand(strand), abundance(abundance), mphfIndex(mphfIndex) , iterationRank(0) {}

    /** kmer value for the node (min of the forward and revcomp value of the bi-directed DB graph). */
    Node_t::Value  kmer;

    /** Strand telling how to interpret the node in the bi-directed DB graph. */
    kmer::Strand strand;

    /** Abundance of the kmer in the initial set of reads. */
    u_int16_t    abundance;

    u_int64_t mphfIndex;
    u_int64_t iterationRank; // used in Simplifications.cpp (see note on tips)

    /** Overload of operator ==  NOTE: by now, it doesn't take care of the strand... */
    bool operator== (const Node_t& other) const  { return kmer == other.kmer; }

    /** Overload of operator !=  NOTE: by now, it doesn't take care of the strand... */
    bool operator!= (const Node_t& other) const  { return kmer != other.kmer; }

    /** Overload of operator <  NOTE: by now, it doesn't take care of the strand... */
    bool operator< (const Node_t& other) const  { return (kmer   < other.kmer); }

    /** Setter for some attributes of the Node object. Usually used by the Graph class for setting
     * a Node object's guts.
     * \param[in] kmer : the kmer value to be set
     * \param[in] strand : strand to be set
     */
    void set (const Node_t::Value& kmer, const kmer::Strand& strand)
    {
        this->kmer      = kmer;
        this->strand    = strand;
        this->mphfIndex = 0;
        this->iterationRank = 0;
    }

    template<typename T>
    const T& getKmer() const;
};

template <> 
template<typename T>
const T& Node_t<tools::math::Integer>::getKmer() const
{
    return this->kmer.get<T>();
}

/********************************************************************************/

/** \brief Specific Node structure representing a branching node in the De Bruijn graph
 *
 * The BranchingNode inherits from the Node structure.
 *
 * Its semantics is to define nodes this way:
 *      indegree!=1 || outdegree!=1
 *
 * Note: with such a definition, nodes that are tips (ie indegree==0 || outdegree=0) are considered
 * as branching. This could be refined if needed (defined this way for historical matters).
 *
 * Its main usage is to be the template specialization type for some Graph class methods.
 */
template<typename Node>
struct BranchingNode_t : Node
{
};

/********************************************************************************/

/** \brief Specific Node structure representing a simple node in the De Bruijn graph
 *
 * The SimpleNode inherits from the Node structure.
 *
 * Its semantics is to define nodes this way:
 *      indegree==1 && outdegree==1
 *
 * Its main usage is to be the template specialization type for some Graph class methods.
 */
template <typename Node>
struct SimpleNode : Node

{
};

/********************************************************************************
                        #######  ######    #####   #######
                        #        #     #  #     #  #
                        #        #     #  #        #
                        #####    #     #  #  ####  #####
                        #        #     #  #     #  #
                        #        #     #  #     #  #
                        #######  ######    #####   #######
********************************************************************************/

/** \brief Definition of an Edge, ie a transition between two nodes in the De Bruijn graph.
 *
 * The Edge structure represents an oriented transition between two nodes. Therefore, it holds:
 *  - the 'from' node
 *  - the 'to' node
 *  - the direction of the transition
 *  - the nucleotide that decorates the transition
 *
 *  The Edge objects are mainly provided by the Graph class, for instance when neighbors of a node
 *  are needed.
 *
 *  The Edge structure may be used as a template specialization type for some Graph class methods.
 */
template <typename Node_t>
struct Edge_t
{
    /** The source node of the edge. */
    Node_t             from;

    /** The target node of the edge. */
    Node_t             to;

    /** The transition nucleotide. */
    kmer::Nucleotide nt;

    /** The direction of the transition. */
    Direction        direction;
    
    /** Overload of operator <.  May not really mean much to compare edges, but is used in Minia's graph simplifications */
    bool operator< (const Edge_t<Node_t>& other) const  { return ((from < other.from) || (from == other.from && to < other.to)); }

    /** Setter for some attributes of the Edge object.
     * \param[in] kmer_from : kmer value of the 'from' Node
     * \param[in] strand_from : strand of the 'from' Node
     * \param[in] kmer_to : kmer value of the 'to' Node
     * \param[in] strand_to : strand of the 'to' Node
     * \param[in] n : the transition nucleotide
     * \param[in] dir : direction of the transition.
     */
    void set (
        const typename Node_t::Value& kmer_from, kmer::Strand strand_from,
        const typename Node_t::Value& kmer_to,   kmer::Strand strand_to,
        kmer::Nucleotide n, Direction dir
    )
    {
        from.set (kmer_from, strand_from);
        to.set   (kmer_to,   strand_to);
        nt = n;
        direction = dir;
    }
};


/********************************************************************************/

/** \brief Specific Edge structure representing a transition between two branching nodes
 *
 * The BranchingEdge inherits from the Edge structure.
 *
 * Note that we don't want here to get the whole simple nodes between the two branching nodes
 * (or all the transition nucleotides required to go from source branching node to target branching node).
 *
 * However, we add a 'distance' attribute that counts the number of transition nucleotides between the two
 * branching nodes.
 */
template<typename Node, typename Edge>
struct BranchingEdge_t : Edge 
{
    void set (
        const typename Node::Value& kmer_from, kmer::Strand strand_from,
        const typename Node::Value& kmer_to,   kmer::Strand strand_to,
        kmer::Nucleotide nt, Direction dir, size_t distance
    )
    {
        Edge::set (kmer_from, strand_from, kmer_to, strand_to, nt, dir);
        this->distance = distance;
    }

    // number of simple nodes between the two branching nodes
    size_t distance;
};

/********************************************************************************
                    ######      #     #######  #     #
                    #     #    # #       #     #     #
                    #     #   #   #      #     #     #
                    ######   #     #     #     #######
                    #        #######     #     #     #
                    #        #     #     #     #     #
                    #        #     #     #     #     #
********************************************************************************/

/** \brief Structure representing a path in the De Bruijn graph.
 *
 * The Path structure provides information on a linear path in the DB graph. It mainly holds
 * the succession of nucleotides that define the path.
 *
 * The path start is defined by an Node object, which defines without ambiguity the first
 * node of the path.
 *
 * It can be used by clients such a contiger tool.
 *
 * Note: by now, the structure is not perfect and some adaptations could come.
 */
template<typename Node>
struct Path_t
{
    /** Constructor (default one)
     * \param[in] n : size of the path. */
    Path_t (size_t n=0) : path(n) {}

    /** Node defining the initial transition of the path. */
    Node start;

    /** Path definition as a succession of nucleotides. */
    std::vector<kmer::Nucleotide> path;

    /** Get the path size
     * \return the size of the path. */
    size_t size() const  { return path.size(); }

    /** Set the size of the path
     * \param[in] n : the size of the path. */
    void resize (size_t n) { path.resize(n); }

    /** Get the ascii value of the ith nucleotide in the path.
     * \param[in] i : index of the nucleotide
     * \return the ascii code for the ith nucleotide. */
    char ascii (size_t i) const { return gatb::core::kmer::ascii((*this)[i]); }

    /** Add a nucleotide to the path.
     * \param[in] nt : nucleotide to be appended to the path. */
    void push_back (kmer::Nucleotide nt)  { path.push_back(nt); }

    /** Clear the path (ie remove all the nucleotides). */
    void clear () { path.clear(); }

    /** Retrieve a reference on the ith nucleotide in the path.
     * \param[in] i : index od the nucleotide to be retrieved.
     * \return a reference on the ith nucleotide */
    kmer::Nucleotide& operator[] (size_t i)        { return path[i]; }

    /** Retrieve a reference on the ith nucleotide in the path.
     * \param[in] i : index od the nucleotide to be retrieved.
     * \return a const reference on the ith nucleotide */
    const kmer::Nucleotide& operator[] (size_t i) const  { return path[i]; }
};

/** Define a comparator for two path. The comparison is a lexicographic comparison on
 * the ascii representation of the path.
 * \param[in] a : path a
 * \param[in] b : path b
 * \return true if path a is less than path b. */
template <typename T>
inline bool operator< (const Path_t<T>& a, const Path_t<T>& b)
{
    size_t N = std::min(a.size(),b.size());
    for (size_t i=0; i<N; i++)
    {
             if (a.ascii(i) < b.ascii(i)) { return true;  }
        else if (a.ascii(i) > b.ascii(i)) { return false; }
    }
    return a.size() < b.size();
}


/** Output stream operator for dumping a Path object as an ascii string
 * holding the nucleotides of the path.
 * \param[in] s : the output stream
 * \param[in] p : the path to be output
 * \return the output stream.
 */
template <typename T>
inline std::ostream& operator<< (std::ostream& s, const Path_t<T>& p)
{
    for (size_t i=0; i<p.size(); i++)  { s << p.ascii(i); }
    return s;
}

/* vector and iterator
 */

/* this was previously a member of GraphTemplate::, and it was named Iterator
 * but I see no reason to, no i'm moving it here and calling it GraphIterator
 * it will help readability and convenience: no more need to pass
 * template parameters <Node,Edge,GraphDataVariant>
 */
template<typename Item>
class GraphIterator : public tools::dp::ISmartIterator<Item>
{
    public:

        /** */
        GraphIterator (tools::dp::ISmartIterator<Item>* ref) : _ref(0) { setRef(ref); }

        /* empty iterator, just for debug purposes, shouldn't be used */
        GraphIterator () : _ref(0) {  std::cout << "empty GraphIterator used (shouldn't happen)" << std::endl; }

        /** */
        virtual ~GraphIterator() { setRef(0); }

        /** */
        GraphIterator (const GraphIterator<Item>& i) : _ref(0) { setRef(i._ref); }

        /** */
        GraphIterator& operator= (const GraphIterator<Item>& i)  {  if (this != &i)   {  setRef (i._ref);  }   return *this;  }

        /** Method that the number of iterated items so far. */
        virtual u_int64_t rank() const  { return _ref->rank(); }

        /** Method that initializes the iteration. */
        void first()  { _ref->first(); }

        /** Method that goes to the next item in the iteration.
         * \return status of the iteration
         */
        void next()  { _ref->next(); }

        /** Method telling whether the iteration is finished or not.
         * \return true if iteration is finished, false otherwise.
         */
        bool isDone()  { return _ref->isDone(); }

        /** Method that returns the current iterated item. Note that the returned type is the template type.
          \return the current item in the iteration.
          */
        Item& item () { return _ref->item(); }

        /** */
        void setItem (Item& i)  {  _ref->setItem (i); }

        /** */
        u_int64_t size () const  { return _ref->size(); }

        /** */
        tools::dp::ISmartIterator<Item>* get()  const { return _ref; }

    private:

        tools::dp::ISmartIterator<Item>* _ref;
        void setRef (tools::dp::ISmartIterator<Item>* ref)  { SP_SETATTR(ref); }
};

/* exactly same comment as GraphIterator above*/
template<typename Item, int NB=16 /* TODO bring it down to 8 whenever possible; or try with a much higher value (e.g. 128) and see if there is actually a performance hit*/>
class GraphVector
{
    public:
        GraphVector () : _size(0)  {} 
        Item& operator[] (size_t idx)  { return (_items[idx]); }
        size_t size()  { return _size; }
        void resize (size_t n)  { _size = n; }

        template<typename Functor>  void iterate (const Functor& f)  { for (size_t i=0; i<_size; i++)  { f(_items[i]); } }

    protected:
        Item   _items[NB];
        size_t _size;
};





/********************************************************************************
                 #####   ######      #     ######   #     #
                #     #  #     #    # #    #     #  #     #
                #        #     #   #   #   #     #  #     #
                #  ####  ######   #     #  ######   #######
                #     #  #   #    #######  #        #     #
                #     #  #    #   #     #  #        #     #
                 #####   #     #  #     #  #        #     #
********************************************************************************/

/** \brief Class representing a De Bruijn graph.
 *
 * This class is the entry point for managing De Bruijn class in gatb-core.
 *
 * Getting a Graph object can be done through :
 *      - creating Graph object (likely from a set of reads).
 *      - loading a Graph object from a file
 *
 * Once a client has a Graph object (with create or load), it is possible to goes
 * through the graph in different ways.
 *
 * The first possibility is to use a Node iterator on the globality of the graph. For
 * instance, all the nodes can be iterated this way, or only branching nodes.
 *
 * The second possibility is to navigate starting from a specific node. For instance,
 * the neighbors of the starting node can be reached.
 *
 * Note: the Graph class doesn't provide means to mark nodes (ie remember which nodes
 * have been visited); this feature could be let to subclasses or other helpers classes.
 *
 * Some utility methods may be useful for debugging (like ascii representation of a node
 * or an edge).
 *
 * The underlying structure of the graph is taken from Minia:
 *      - a Bloom filter
 *      - a set of false positives
 *
 * Once a graph is built (from a set of reads), it is saved in a file (likely HDF5 format).
 * It is so possible to get a Graph object by loading the file instead of re-build it.
 *
 * Note: branching nodes are computed during the graph building; they are also saved in the graph
 * output file.
 */
template <typename Node, typename Edge, typename GraphDataVariant>
class GraphTemplate
{
public:

    /********************************************************************************/
    /*                            STATIC METHODS   (create/load)                    */
    /********************************************************************************/

    /** Build an empty graph.
     * \param[in] kmerSize: kmer size
     * \return the created graph.
     */
    static GraphTemplate  create (size_t kmerSize)  {  return  GraphTemplate(kmerSize);  }

    /** Build a graph from a given bank.
     * \param[in] bank : bank to get the reads from
     * \param[in] fmt : printf-like format for the command line string
     * \return the created graph.
     */
    static GraphTemplate  create (bank::IBank* bank, const char* fmt, ...);

    /** Build a graph from user options.
     * \param[in] fmt: printf-like format
     * \return the created graph.
     */
    static GraphTemplate  create (const char* fmt, ...);

    /** Build a graph from scratch.
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static GraphTemplate  create (tools::misc::IProperties* options)  {  return  GraphTemplate (options);  }

    /** Load a graph from some URI.
     * \param[in] uri : the uri to get the graph from
     * \return the loaded graph.
     */
    static GraphTemplate  load (const std::string& uri)  {  return  GraphTemplate (uri);  }

    /** Get a parser object that knows the user options for building a graph.
     * \return the options parser object.
     */
    static tools::misc::IOptionsParser* getOptionsParser (bool includeMandatory=true);

    /********************************************************************************/
    /*                               CONSTRUCTORS                                   */
    /********************************************************************************/

    /* Default Constructor.*/
    GraphTemplate ();

    /* Copy Constructor.*/
    GraphTemplate (const GraphTemplate& graph);

    /* Destructor. */
    ~GraphTemplate ();

    /** Affectation overload. */
    GraphTemplate& operator= (const GraphTemplate& graph);

    /**********************************************************************/
    /*                     GLOBAL ITERATOR METHODS                        */
    /**********************************************************************/

    /** Creates an iterator over nodes of the graph.
     *      this used to be a templated method but I'm now untemplating it, because of nested templates specialization
     *      so call iteratorBranching if you want an iterator over BranchingNode's
     * \return the nodes iterator. */

    inline GraphIterator<Node> iterator () const  {  return getNodes ();           }
    inline GraphIterator<BranchingNode_t<Node> > iteratorBranching () const  {  return getBranchingNodes ();           }
    GraphIterator<Node> iteratorCachedNodes () const;


    /**********************************************************************/
    /*                     ALL NEIGHBORS METHODS                          */
    /**********************************************************************/

    /** Returns a vector of neighbors of the provided node.
     * \param[in] node : the node whose neighbors are wanted
     * \param[in] direction : the direction of the neighbors. If not set, out and in neighbors are computed.
     * \return a vector of the node neighbors (may be empty). 
     * Warning: be sure to check if edge.from (or node.from) is actually your input node, or its reverse complement.
     */
    inline GraphVector<Node> neighbors    ( Node& node, Direction dir=DIR_END) const  {  return getNodes(node, dir);           }
    
    /** Returns a vector of neighbors of the provided kmer. It has to be understood as the following:
     *  - a node N is built with the kmer, with the strand FORWARD
     *  - a call to 'neighbors<T> (N,          DIR_OUTGOING)' is done; we get v1
     *  - a call to 'neighbors<T> (reverse(N), DIR_OUTGOING)' is done; we get v2
     *  - the result is the concatenation of v1 and v2
     *  \param[in] kmer : the kmer whose neighbors are wanted.
     *  \return a vector of the neighbors (may be empty).
     */
    inline GraphVector<Node> neighbors    ( const typename Node::Value& kmer) const          {  return getNodeValues (kmer);           }
    
    inline Node* neighborsDummy      ( Node& node, Direction dir=DIR_END) const  {   return NULL;           }


    /* neighbors used to be templated.. not anymore, trying to avoid nested template specialization; 
     * so call neighbors for getting nodes, or neighborsEdge for getting edges */
    inline GraphVector<Edge> neighborsEdge    ( Node& node, Direction dir=DIR_END) const  {   return getEdges(node, dir);           }
    inline Edge* neighborsDummyEdge      ( Node& node, Direction dir=DIR_END) const  {   return NULL;           }
    inline GraphVector<Edge> neighborsEdge    ( const typename Node::Value& kmer) const          {  return getEdgeValues (kmer);           }


    inline GraphVector<BranchingNode_t<Node> > neighborsBranching(Node& node, Direction direction) const
    { return getBranchingNodeNeighbors (node, direction);  }

    inline GraphVector<BranchingNode_t<Node> > neighborsBranching(const typename Node::Value& kmer) const
    { return getBranchingNodeValues (kmer);  }

     inline GraphVector<BranchingEdge_t<Node,Edge> > neighborsBranchingEdge (Node& node, Direction direction) const
    { return getBranchingEdgeNeighbors (node, direction);  }

    inline GraphVector<BranchingEdge_t<Node,Edge> > neighborsBranchingEdge (const typename Node::Value& kmer) const
    { return getBranchingEdgeValues (kmer);  }



    /** Shortcut for 'neighbors' method with direction==DIR_OUTCOMING.
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    inline GraphVector<Node> successors   ( Node& node) const                 {  return getNodes(node, DIR_OUTCOMING); }

    /** Shortcut for 'neighbors' method with direction==DIR_INCOMING.
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    inline GraphVector<Node> predecessors ( Node& node) const                 {  return getNodes(node, DIR_INCOMING);  }
    
    inline GraphVector<Edge> successorsEdge   ( Node& node) const                 {  return getEdges(node, DIR_OUTCOMING); }
    inline GraphVector<Edge> predecessorsEdge ( Node& node) const                 {  return getEdges(node, DIR_INCOMING);  }

    inline GraphVector<BranchingNode_t<Node> > successorsBranching (Node& node) const
        { return getBranchingNodeNeighbors (node, DIR_OUTCOMING);  }

    inline GraphVector<BranchingNode_t<Node> > predecessorsBranching (Node& node) const
        { return getBranchingNodeNeighbors (node, DIR_INCOMING);  }

    inline GraphVector<BranchingEdge_t<Node,Edge> > successorsBranchingEdge (Node& node) const
        { return getBranchingEdgeNeighbors (node, DIR_OUTCOMING);  }

    inline GraphVector<BranchingEdge_t<Node,Edge> > predecessorsBranchingEdge (Node& node) const
        { return getBranchingEdgeNeighbors (node, DIR_INCOMING);  }


    /** Returns a set of neighbors for each node iterated with the provided two iterators
     * \param[in] first : beginning of the iteration
     * \param[in] last : end of the iteration
     * \return all the neighbors computed for each iterated node. */
    /* used to be a template, but I can't specialize it without specializing the whole graph templated class now. so got rid of templates for now*/
    //template <typename T, typename IteratorInput>
    //std::set<T> neighbors (IteratorInput first, IteratorInput last) const;
    std::set<BranchingNode_t<Node> > neighbors (typename std::set<BranchingNode_t<Node> >::iterator first, typename std::set<BranchingNode_t<Node> >::iterator last) const;

    /** Returns the successors of two nodes, ie with the same transition nucleotide from both nodes.
     * \param[in] node1 : first node
     * \param[in] node2 : sedond node
     * \return the vector of pairs of items as successors
     */
    inline GraphVector<std::pair<Edge,Edge> > successorsEdge (const Node& node1, const Node& node2) const
    { return getEdgesCouple (node1, node2, DIR_OUTCOMING); }

    inline GraphVector<std::pair<Node,Node> > successors (const Node& node1, const Node& node2) const
    { return getNodesCouple (node1, node2, DIR_OUTCOMING); }

    /** Returns the predecessors of two nodes, ie with the same transition nucleotide from both nodes.
     * \param[in] node1 : first node
     * \param[in] node2 : sedond node
     * \return the vector of pairs of items as predecessors
     */
    inline GraphVector<std::pair<Edge,Edge> > predecessorsEdge (const Node& node1, const Node& node2) const
    { return getEdgesCouple (node1, node2, DIR_INCOMING); }

    inline GraphVector<std::pair<Node,Node> > predecessors (const Node& node1, const Node& node2) const
    { return getNodesCouple (node1, node2, DIR_INCOMING); }


    /**********************************************************************/
    /*                     ONE NEIGHBOR METHODS                           */
    /**********************************************************************/

    /** Return a specific neighbor from a given node. The neighbor is defined by a direction and the transition
     * nucleotide.
     * IMPORTANT: this method will not check that the neighbor node belongs to the graph: it merely
     * computes the next kmer but doesn't check the Bloom filter. It is supposed that the client has already asked
     * for the neighbors and so knows the valid transitions.
     * \param[in] source : the source neighbor
     * \param[in] dir : the direction of the transition
     * \param[in] nt : the nucleotide of the transition
     * \return the neighbor object.
     */
    inline Node neighbor ( Node& source, Direction dir, kmer::Nucleotide nt) const
    {  bool exists=true; return getNode (source, dir, nt, exists);  }

    /** Return a specific neighbor from a given node. The neighbor is defined by a direction and the transition
     * nucleotide.
     * IMPORTANT: this method will check that the neighbor node belongs to the graph. If the neighbor is not in the
     * graph, the 'exists' parameter is set to false, true otherwise.
     * \param[in] source : the source neighbor
     * \param[in] dir : the direction of the transition
     * \param[in] nt : the nucleotide of the transition
     * \param[out] exists : yes means that the neighbor is in the graph, false otherwise
     * \return the neighbor object.
     */
    inline Node neighbor ( Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
    {  return getNode (source, dir, nt, exists);  }

    /** Shortcut for neighbor with dir==DIR_OUTCOMING. */
    inline Node successor ( Node& source, kmer::Nucleotide nt) const
    {  bool exists=true; return getNode (source, DIR_OUTCOMING, nt, exists);  }

    inline Node successor ( Node& source, kmer::Nucleotide nt, bool& exists) const
    {  return getNode (source, DIR_OUTCOMING, nt, exists);  }

   
    /** Shortcut for neighbor with dir==DIR_INCOMING. */
    inline Node predecessor ( Node& source, kmer::Nucleotide nt) const
    {  bool exists=true; return getNode (source, DIR_INCOMING, nt, exists);  }

    inline Node predecessor ( Node& source, kmer::Nucleotide nt, bool& exists) const
    {  return getNode (source, DIR_INCOMING, nt, exists);  }


    /**********************************************************************/
    /*                      MISC NEIGHBORS METHODS                        */
    /**********************************************************************/

    /** Get the incoming degree of the node.
     * \param[in] node : the node
     * \return the indegree of the node. */
    size_t indegree  (Node& node) const;

    /** Get the outcoming degree of the node.
     * \param[in] node : the node
     * \return the outdegree of the node. */
    size_t outdegree (Node& node) const;

    /** Get the degree of the node (either incoming or outcoming).
     * \param[in] node : the node
     * \param[in] dir : direction of the degree
     * \return the degree of the node. */
    size_t degree    (Node& node, Direction dir) const;
    
    /* get both the in and out degree of a node, in a single call
     * takes advantage of adjacency
     */
    void degree    (Node& node, size_t& in, size_t &out) const;
    
    /**********************************************************************/
    /*                      SIMPLIFICATION METHODS                        */
    /**********************************************************************/

    /* perform tip removal, bulge removal and EC removal, as in Minia */
    void simplify(unsigned int nbCores = 1, bool verbose=true);

    /**********************************************************************/
    /*                         SIMPLE PATH METHODS                        */
    /**********************************************************************/

    /** Simple paths traversal
     *  invariant: the input kmer has no in-branching.
     * \returns
     *       1 if a good extension is found
     *       0 if a deadend was reached
     *      -1 if out-branching was detected
     *      -2 if no out-branching but next kmer has in-branching
     */
    int simplePathAvance (Node& node, Direction dir, Edge& output) const;

    /** */
    int simplePathAvance (Node& node, Direction dir) const;

    /** */
    int simplePathAvance (Node& node, Direction dir, kmer::Nucleotide& nt) const;

    /** */
    GraphIterator<Node> simplePath     (Node& node, Direction dir) const  { return getSimpleNodeIterator(node, dir); }
    GraphIterator<Edge> simplePathEdge (Node& node, Direction dir) const  { return getSimpleEdgeIterator(node, dir); }
    
    // high-level functions that used to be in Simplifications.cpp
    // convention: unitigXXX works on the unitigs as computed as bcalm. never leaves that unitig
    //             simplepathXXX may traverse multiple unitigs
    Node             unitigLastNode          (Node& node, Direction dir) const;
    Node         simplePathLastNode          (Node& node, Direction dir) const;
    unsigned int     unitigLength            (Node& node, Direction dir) const;
    unsigned int simplePathLength            (Node& node, Direction dir) const;
    double           unitigMeanAbundance     (Node& node) const;
    double       simplePathMeanAbundance     (Node& node, Direction dir) const;
    void             unitigDelete          (Node& node, Direction dir, NodesDeleter<Node,Edge, GraphTemplate<Node, Edge, GraphDataVariant> >& nodesDeleter) ;
    void             unitigDelete          (Node& node) ;
    void         simplePathDelete          (Node& node, Direction dir, NodesDeleter<Node,Edge, GraphTemplate<Node, Edge, GraphDataVariant> >& nodesDeleter) ;
    std::string  unitigSequence (Node& node, bool& isolatedLeft, bool& isolatedRight) const;
    std::string  unitigPathSequence (Node& node, bool& isolatedLeft, bool& isolatedRight) const;
    void         unitigMark            (Node& node) ; // used to flag simple path as traversed, in minia
    bool         unitigIsMarked        (Node& node) ;
    
    std::string simplePathBothDirections(const Node& node, bool& isolatedLeft, bool& isolatedRight, bool dummy, float& coverage);
    // aux function, not meant to be called from outside, but maybe it could.
    void simplePathLongest_avance(const Node& node, Direction dir, int& seqLength, int& endDegree, bool markDuringTraversal, float& coverage, std::string* seq = nullptr, std::vector<Node> *unitigNodes = nullptr); 
    
    void debugPrintAllUnitigs() const;

    /**********************************************************************/
    /*                         NODE METHODS                               */
    /**********************************************************************/

    /** Tells whether or not a node belongs to the graph.
     * \param[in] item : the node
     * \return true if the node belongs to the graph, false otherwise. */
    bool contains (const Node& item) const;
                                  
    template<size_t span> 
    bool contains (const typename gatb::core::kmer::impl::Kmer<span>::Type& item) const;

    /** Get the ascii string for the node, according to its strand.
     * \param[in] node: the node to get the string from
     * \return the string representation for the provided node. */
    std::string toString (const Node& node) const;

    /** Tells whether the provided node is branching or not.
     * \param[in] node : the node to be asked
     * \return true if the node is branching, false otherwise. */
    bool isBranching (Node& node) const;

    /** Build a fake node (ie. not necessarily in the De Bruijn graph). Mainly for test purpose.
     * \param[in] data : a string like structure for the sequence from which the kmer of the node is extracted
     * \param[in] offset : starting offset in the data
     * \return the fake node. */
    Node buildNode (const tools::misc::Data& data, size_t offset=0) const;

    /** Build a fake node (ie. not necessarily in the De Bruijn graph). Mainly for test purpose.
     * \param[in] sequence : a sequence of nucleotides in ASCII format
     * \return the fake node. */
    Node buildNode (const char* sequence) const;

    /** Return the reverse complement node of the provided one.
     * param[in] node : the node to be reverted
     * \return the reverted node.  */
    Node reverse (const Node& node) const;

    /** Return the reverse complement node of the provided one.
     * param[in] node : the node to be reverted
     * \return the reverted node.  */
    BranchingNode_t<Node> reverse (const BranchingNode_t<Node>& node) const;

    /** Mutation of a node. */
    GraphVector<Node> mutate (const Node& node, size_t idx, int mode=0) const;

    /** Return a nucleotide at position 'idx' of a given node 
     * \param[in] node : the node we want to extract a nucleotide from
     * \param[in] idx : the position of the nucleotide to be extracted
     * \return the wanted nucleotide. */
    kmer::Nucleotide getNT (const Node& node, size_t idx) const;

    /** Return the abundance of a node by querying the perfect hash function 
     * \param[in] node : the node
     * \return the abundance */
    int queryAbundance (Node& node) const;

    /** Return the state of a node by querying the perfect hash function. A node state is either normal, marked, or deleted.
     * \param[in] node : the node or a node index (unsigned long) from the MPHF
     * \return the abundance */
    int queryNodeState (Node& node) const;
    void setNodeState (Node& node, int state) const;
    void resetNodeState () const ;
    void disableNodeState () const ; // see Graph.cpp for explanation

    // deleted nodes, related to NodeState above
    void deleteNode (Node& node) const;
    void deleteNodesByIndex(std::vector<bool> &bitmap, int nbCores = 1, gatb::core::system::ISynchronizer* synchro=NULL) const;
    bool isNodeDeleted(Node& node) const;

    // a direct query to the MPHF data strcuture
    /* returns an index between 0 and (number of nodes - 1)
     * the function is a bijection between the nodes and the indices
     *
     * NOTE: if you query this function for a node that doesn't beling to the graph, 
     * it may still return an index (that will be in collision with
     * the index of another node in the graph), or it may return ULLONG_MAX, 
     * the latter indicates that for sure, the node isn't in the graph
     */
    unsigned long nodeMPHFIndex(Node& node) const;


    unsigned long nodeMPHFIndexDummy(Node& node) const; // debug function, for profiling only


    /**********************************************************************/
    /*                         EDGE METHODS                               */
    /**********************************************************************/

    /** Get the ascii string for the edge
     * \param[in] edge : the edge to get the string from
     * \return the string representation for the provided edge . */
    std::string toString (const Edge& edge) const;

    /** Get the ascii string for the branching edge
     * \param[in] edge : the edge to get the string from
     * \return the string representation for the provided edge . */
    std::string toString (const BranchingEdge_t<Node, Edge>& edge) const;

    /** Tells whether the provided edge is simple: outdegree(from)==1 and indegree(to)==1
     * \param[in] edge : the edge to be asked
     * \return true if the edge is simple, false otherwise. */
    bool isSimple (Edge& edge) const;

    /**********************************************************************/
    /*                         MISC METHODS                               */
    /**********************************************************************/

    /** Return the name of the graph.
     * \return the name. */
    std::string getName() const { return _name; }

    /** Get the size of the kmers.
     * \return the kmer size. */
    size_t getKmerSize() const { return _kmerSize; }

    /** Get information about the graph (gathered during its creation).
     * \return a property object holding graph information. */
    tools::misc::IProperties& getInfo () const { return (tools::misc::IProperties&)_info; }

    /** Remove physically a graph. */
    void remove ();

    /** Reverse an edge.
     * param[in] edge: the edge to be reverted
     * \return the reverted edge. */
    Edge reverse (const Edge& edge) const;


    /** cache adjacency information from the Bloom filter to an array, 8 bits per node, for faster traversal queries*/
    void precomputeAdjacency(unsigned int nbCores = 1, bool verbose = true);
    unsigned int nt2bit[256]; 
    bool debugCompareNeighborhoods(Node& node, Direction dir, std::string prefix) const; // debug


    /* functions to precompute non-simple nodes */
    void cacheNonSimpleNode(const Node& node) const;
    void cacheNonSimpleNodeDelete(const Node& node) const;
    void cacheNonSimpleNodes(unsigned int nbCores, bool verbose) ;

    /**********************************************************************/
    /*                         DEBUG METHODS                              */
    /**********************************************************************/
    /** */
    std::string debugString (const Node& node, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const;

    /** */
    std::string debugString (const Edge& edge, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const;


    /**********************************************************************/
    /*                         TYPES                                      */
    /**********************************************************************/
    enum StateMask
    {
        STATE_INIT_DONE           = (1<<0),
        STATE_CONFIGURATION_DONE  = (1<<1),
        STATE_SORTING_COUNT_DONE  = (1<<2),
        STATE_BLOOM_DONE          = (1<<3),
        STATE_DEBLOOM_DONE        = (1<<4),
        STATE_BRANCHING_DONE      = (1<<5),
        STATE_MPHF_DONE           = (1<<6),
        STATE_ADJACENCY_DONE      = (1<<7),
        STATE_NONSIMPLE_CACHE     = (1<<8)
    };
    typedef u_int64_t State; /* this is a global graph state, not to be confused of the state of a node (deleted or not) */
    State getState () const { return _state; }
    bool  checkState (StateMask mask) const { return (_state & (State)mask)==(State)mask; }
    State setState   (StateMask mask) { _state |=  mask; return _state; }
    State unsetState (StateMask mask) { _state &= ~mask; return _state; }

    /** Constructor for empty graph.*/
    GraphTemplate (size_t kmerSize);

    /** Constructor. Use for GraphTemplate creation (ie. DSK + debloom) and filesystem save. */
    GraphTemplate (bank::IBank* bank, tools::misc::IProperties* params);

    /** Constructor. Use for GraphTemplate creation (ie. DSK + debloom) and filesystem save. */
    GraphTemplate (tools::misc::IProperties* params);

    /** Constructor. Use for reading from filesystem. */
    GraphTemplate (const std::string& uri);

public: // was private: before, but had many compilation errors during the change from Graph to GraphTemplate. took the easy route, set it to "public:", it solved everything.
    
    /** Kind of storage for the graph. */
    tools::storage::impl::StorageMode_e _storageMode;

    /** Default storage kind. */
    static const tools::storage::impl::StorageMode_e PRODUCT_MODE_DEFAULT = tools::storage::impl::STORAGE_HDF5;

    /** Storage. */
    tools::storage::impl::Storage* _storage;
    void setStorage (tools::storage::impl::Storage* storage)  { SP_SETATTR(storage); }
    tools::storage::impl::Storage& getStorage()                           { return (*_storage); }
    tools::storage::impl::Group&   getGroup  (const std::string name="")  { return getStorage() (name); }

    /** Defined as a void* for hiding implementation in cpp file. */
    void* _variant;

    /** kmer size of the graph */
    size_t _kmerSize;

    /** Creation information. */
    tools::misc::impl::Properties _info;

    /** */
    std::string _name;

    State _state;

    /** */
    tools::misc::BloomKind       _bloomKind;
    tools::misc::DebloomKind     _debloomKind;
    tools::misc::DebloomImpl     _debloomImpl;
    tools::misc::BranchingKind   _branchingKind;
   
    /** */
    GraphIterator<Node> getNodes () const;
    
    unsigned char countNeighbors (Node&, Direction) const; // simple and much faster version of getNodes, for degree(), outdegree(), indegree() queries
    void countNeighbors (Node&, size_t&, size_t&) const;  // compute in and out degree at the same time

    /** */
    GraphIterator<BranchingNode_t<Node> > getBranchingNodes () const;

    /** */
    GraphIterator<Node> getSimpleNodeIterator (Node& node, Direction dir) const;

    /** */
    GraphIterator<Edge> getSimpleEdgeIterator (Node& node, Direction dir) const;

    /** */
    GraphVector<Edge> getEdges (Node source, Direction direction) const;

    /** */
    GraphVector<std::pair<Node,Node> > getNodesCouple (const Node& node1, const Node& node2, Direction direction) const;

    /** */
    GraphVector<std::pair<Edge,Edge> > getEdgesCouple (const Node& node1, const Node& node2, Direction direction) const;

    /** */
    GraphVector<Node> getNodes (Node &source, Direction direction)  const;

    /** */
    GraphVector<BranchingNode_t<Node> > getBranchingNodeNeighbors (Node& source, Direction direction) const;

    /** */
    GraphVector<BranchingEdge_t<Node,Edge> > getBranchingEdgeNeighbors (Node& source, Direction direction) const;

    /** */
    GraphVector<Edge> getEdgeValues (const typename Node::Value& kmer) const;

    /** */
    GraphVector<Node> getNodeValues (const typename Node::Value& kmer) const;

    /** */
    GraphVector<BranchingEdge_t<Node,Edge> > getBranchingEdgeValues (const typename Node::Value& kmer) const;

    /** */
    GraphVector<BranchingNode_t<Node> > getBranchingNodeValues (const typename Node::Value& kmer) const;

    /** */
    Node getNode (Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const;

    /* set the graph variant */
    void setVariant (void* data, size_t kmerSize, size_t integerPrecision=0);

    /** Friends. */
    template<typename, typename, typename> friend struct build_visitor_solid ; // i don't know why this template<typename, typename> trick works, but it does
    template<typename, typename, typename> friend struct build_visitor_postsolid ;
    template<typename, typename, typename> friend struct configure_visitor;

    // a late addition, because GraphUnitig wants to call it too
    static void executeAlgorithm (gatb::core::tools::misc::impl::Algorithm& algorithm, gatb::core::tools::storage::impl::Storage* storage, gatb::core::tools::misc::IProperties* props, gatb::core::tools::misc::IProperties& info);
};


/********************************************************************************
                        #     #  ###   #####    #####
                        ##   ##   #   #     #  #     #
                        # # # #   #   #        #
                        #  #  #   #    #####   #
                        #     #   #         #  #
                        #     #   #   #     #  #     #
                        #     #  ###   #####    #####
********************************************************************************/

template<class Type, class Listener>
class ProgressGraphIteratorTemplate : public tools::dp::impl::SubjectIterator<Type>
{
public:
    ProgressGraphIteratorTemplate (const GraphIterator<Type>& items, const char* msg = "compute", bool verbose=true, size_t divide=100)
        : tools::dp::impl::SubjectIterator<Type> (items.get(), items.size()/divide), _size(items.size()) 
        { 
            if (verbose)
                this->addObserver( new Listener(items.size(), msg) );
        }
    

    u_int64_t size () const  { return _size; }

private:
    u_int64_t _size;
};

// I tried to overload std::hash in Integer but it didn't work, despite the fact that it should have!
template<typename T>
class NodeHasher
{
    public:
        size_t operator()(const T & k) const 
        {
            return hash1(k, 0);
        }
};

/********************************************************************************/

/* We define a structure that holds all the necessary stuff for implementing the graph API.
 *  Here, the structure is templated by the span (ie. the kmer max size).
 *
 *  This structure is the basis for defining a boost::variant with all required span
 *  template instantiations.
 */
template<size_t span>
struct GraphData
{
    /** Shortcuts. */
    typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename gatb::core::kmer::impl::Kmer<span>::Type           Type;
    typedef typename gatb::core::kmer::impl::Kmer<span>::Count          Count;
    typedef typename gatb::core::kmer::impl::MPHFAlgorithm<span>::AbundanceMap   AbundanceMap;
    typedef typename gatb::core::kmer::impl::MPHFAlgorithm<span>::NodeStateMap   NodeStateMap;
    typedef typename gatb::core::kmer::impl::MPHFAlgorithm<span>::AdjacencyMap   AdjacencyMap;
    typedef typename std::unordered_map<Type, std::pair<char,std::string>, NodeHasher<Type> > NodeCacheMap; // rudimentary for now

    /** Constructor. */
    GraphData () : _model(0), _solid(0), _container(0), _branching(0), _abundance(0), _nodestate(0), _adjacency(0), _nodecache(0) {}

    /** Destructor. */
    ~GraphData ()
    {
        setModel     (0);
        setSolid     (0);
        setContainer (0);
        setBranching (0);
        setAbundance (0);
        setNodeState (0);
        setAdjacency (0);
        setNodeCache (0);
    }

    /** Constructor (copy). */
    GraphData (const GraphData& d) : _model(0), _solid(0), _container(0), _branching(0), _abundance(0), _nodestate(0), _adjacency(0), _nodecache(0)
    {
        setModel     (d._model);
        setSolid     (d._solid);
        setContainer (d._container);
        setBranching (d._branching);
        setAbundance (d._abundance);
        setNodeState (d._nodestate);
        setAdjacency (d._adjacency);
        setNodeCache (d._nodecache);
    }

    /** Assignment operator. */
    GraphData& operator= (const GraphData& d)
    {
        if (this != &d)
        {
            setModel     (d._model);
            setSolid     (d._solid);
            setContainer (d._container);
            setBranching (d._branching);
            setAbundance (d._abundance);
            setNodeState (d._nodestate);
            setAdjacency (d._adjacency);
            setNodeCache (d._nodecache);
        }
        return *this;
    }

    /** Required attributes. */
    Model*                _model;
    tools::storage::impl::Partition<Count>*   _solid;
    IContainerNode<Type>*                     _container;
    tools::collections::Collection<Count>*    _branching;
    AbundanceMap*         _abundance;
    NodeStateMap*         _nodestate;
    AdjacencyMap*         _adjacency;
    NodeCacheMap*         _nodecache; // so, nodecache also records branching node, but also more stuff. i'm keeping _branching for historical reasons.

    /** Setters. */
    void setModel       (Model*                                       model)      { SP_SETATTR (model);     }
    void setSolid       (tools::storage::impl::Partition<Count>*      solid)      { SP_SETATTR (solid);     }
    void setContainer   (IContainerNode<Type>*                    container)  { SP_SETATTR (container); }
    void setBranching   (tools::collections::Collection<Count>*   branching)  { SP_SETATTR (branching); }
    void setAbundance   (AbundanceMap*          abundance)  { SP_SETATTR (abundance); }
    void setNodeState   (NodeStateMap*          nodestate)  { SP_SETATTR (nodestate); }
    void setAdjacency   (AdjacencyMap*          adjacency)  { SP_SETATTR (adjacency); }
    void setNodeCache   (NodeCacheMap*          nodecache)  { _nodecache = nodecache; /* would like to do "SP_SETATTR (nodecache)" but nodecache is an unordered_map, not some type that derives from a smartpointer. so one day, address this. I'm not sure if it's important though. Anyway I'm phasing out NodeCache in favor of GraphUnitigs. */; }

    /** Shortcut. */
    bool contains (const Type& item)  const  {  

        bool res = _container->contains (item);

        if (!res)
            return false;

        /* check if kmer is deleted*/
        // this is duplicated code from queryNodeState.
        // NOTE: this does a MPHF query for each bloom contains that answer true. costly!
        if (_nodestate != NULL)
        {
            unsigned long hashIndex = ((_nodestate))->getCode(item);
			if(hashIndex == ULLONG_MAX) return false;
            unsigned char value = ((_nodestate))->at(hashIndex / 2);
            if ((hashIndex % 2) == 1)
                value >>= 4;
            value &= 0xF;
            if (((value >> 1) & 1) == 1) 
                return false;
        }

        return true;
    }
};

/* This definition is the basis for having a "generic" Graph class, ie. not relying on a template
 * parameter.
 *
 * This is done through a boost::variant; actually, we use a limited number of variant, corresponding
 * to some maximum kmer sizes.
 */
template<typename T>  struct ToGraphDataVariant  {  typedef GraphData<T::value> type;  };
/********************************************************************************/


/* typedef for compatibility with all existing GATB tools */

typedef Node_t<> Node; // default Node type: tools::math::Integer (which is a boost::variant)
typedef Edge_t<Node> Edge;
typedef Path_t<Node> Path; 
typedef BranchingNode_t<Node> BranchingNode; 
typedef BranchingEdge_t<Node,Edge> BranchingEdge; 
typedef boost::make_variant_over<boost::mpl::transform<gatb::core::tools::math::IntegerList, ToGraphDataVariant<boost::mpl::_> >::type >::type GraphDataVariant;
typedef GraphTemplate<Node, Edge, GraphDataVariant> Graph; // define classical GATB Graph


/* rationale: Node is when you have no idea what the kmer size is going to be. NodeFast is when you do. and it's faster */

template <size_t span>
using NodeFast = Node_t<typename gatb::core::kmer::impl::Kmer<span>::Type >;
template <size_t span>
using EdgeFast = Edge_t<NodeFast<span> >;
template <size_t span>
using GraphDataVariantFast = boost::variant<GraphData<span> >; 


template <typename Type, class Listener>
class ProgressGraphIterator : public ProgressGraphIteratorTemplate<Type, Listener> {
    public:
    ProgressGraphIterator (const GraphIterator<Type>& items, const char* msg = "compute", size_t divide=100) : 
        ProgressGraphIteratorTemplate<Type, Listener> (items,msg,divide)
    {}
};
    

template<typename Node, typename Edge, typename GraphDataVariant>
struct configure_visitor : public boost::static_visitor<>    {

    const GraphTemplate<Node, Edge, GraphDataVariant>& graph;
    tools::storage::impl::Storage&     storage;

    configure_visitor (const GraphTemplate<Node, Edge, GraphDataVariant>& graph, tools::storage::impl::Storage& storage)  : graph(graph), storage(storage) {}

    template<size_t span>  void operator() (GraphData<span>& data) const;
};

template<typename Node, typename Edge, typename GraphDataVariant>
struct build_visitor_solid : public boost::static_visitor<>    {

    GraphTemplate<Node, Edge, GraphDataVariant>& graph; 
    bank::IBank* bank; 
    tools::misc::IProperties* props;

    build_visitor_solid (GraphTemplate<Node, Edge, GraphDataVariant>& aGraph, bank::IBank* aBank, tools::misc::IProperties* aProps)  : graph(aGraph), bank(aBank), props(aProps) {}

    template<size_t span>  void operator() (GraphData<span>& data) const;
};

/* now build the rest of the graph */
template<typename Node, typename Edge, typename GraphDataVariant>
struct build_visitor_postsolid : public boost::static_visitor<>    {

    GraphTemplate<Node, Edge, GraphDataVariant>& graph; tools::misc::IProperties* props; 

    build_visitor_postsolid (GraphTemplate<Node, Edge, GraphDataVariant>& aGraph, tools::misc::IProperties* aProps)  : graph(aGraph), props(aProps) {}

    template<size_t span>  void operator() (GraphData<span>& data) const;
};


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
