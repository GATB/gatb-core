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

#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/tools/math/Integer.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <vector>
#include <set>

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
    DIR_INCOMING,
    DIR_END
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
 *  the initial set of reads. Note that this attribute has a value only for so-called 'solid kmers',
 *  generally computed by the DSK tool.
 */
struct Node
{
    /** Type for a kmer value. */
    typedef tools::math::Integer Value;

    /** Default constructor. */
    Node () : abundance(0), strand(kmer::STRAND_FORWARD) {}

    /** Constructor.
     * \param[in] kmer : kmer value. By default, it is the minimum value of the forward and revcomp value.
     * \param[in] strand : strand telling how to interpret the kmer value. Default value is forward by convention.
     * \param[in] abundance : abundance of the kmer. Default value is 0 if not set.
     */
    Node (const Node::Value& kmer, kmer::Strand strand=kmer::STRAND_FORWARD, u_int16_t abundance=0)
        : kmer(kmer), strand(strand), abundance(abundance) {}

    /** kmer value for the node (min of the forward and revcomp value of the bi-directed DB graph). */
    Node::Value  kmer;

    /** Strand telling how to interpret the node in the bi-directed DB graph. */
    kmer::Strand strand;

    /** Abundance of the kmer in the initial set of reads. */
    u_int16_t    abundance;

    /** Overload of operator ==  NOTE: by now, it doesn't take care of the strand... */
    bool operator== (const Node& other) const  { return kmer == other.kmer; }

    /** Overload of operator !=  NOTE: by now, it doesn't take care of the strand... */
    bool operator!= (const Node& other) const  { return kmer != other.kmer; }

    /** Overload of operator <  NOTE: by now, it doesn't take care of the strand... */
    bool operator< (const Node& other) const  { return (kmer   < other.kmer); }

    /** Setter for some attributes of the Node object. Usually used by the Graph class for setting
     * a Node object's guts.
     * \param[in] kmer : the kmer value to be set
     * \param[in] strand : strand to be set
     */
    void set (const Node::Value& kmer, const kmer::Strand& strand)
    {
        this->kmer      = kmer;
        this->strand    = strand;
    }
};

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
struct BranchingNode : Node
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
struct Edge
{
    /** The source node of the edge. */
    Node             from;

    /** The target node of the edge. */
    Node             to;

    /** The transition nucleotide. */
    kmer::Nucleotide nt;

    /** The direction of the transition. */
    Direction        direction;

    /** Setter for some attributes of the Edge object.
     * \param[in] kmer_from : kmer value of the 'from' Node
     * \param[in] strand_from : strand of the 'from' Node
     * \param[in] kmer_to : kmer value of the 'to' Node
     * \param[in] strand_to : strand of the 'to' Node
     * \param[in] n : the transition nucleotide
     * \param[in] dir : direction of the transition.
     */
    void set (
        const Node::Value& kmer_from, kmer::Strand strand_from,
        const Node::Value& kmer_to,   kmer::Strand strand_to,
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
struct BranchingEdge : Edge
{
    void set (
        const Node::Value& kmer_from, kmer::Strand strand_from,
        const Node::Value& kmer_to,   kmer::Strand strand_to,
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
struct Path
{
    /** Constructor (default one)
     * \param[in] n : size of the path. */
    Path (size_t n=0) : path(n) {}

    /** Edge defining the initial transition of the path. */
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
inline bool operator< (const Path& a, const Path& b)
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
inline std::ostream& operator<< (std::ostream& s, const Path& p)
{
    for (size_t i=0; i<p.size(); i++)  { s << p.ascii(i); }
    return s;
}

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
class Graph
{
public:

    /********************************************************************************/
    /*                            STATIC METHODS   (create/load)                    */
    /********************************************************************************/

    /** Build an empty graph.
     * \param[in] kmerSize: kmer size
     * \return the created graph.
     */
    static Graph  create (size_t kmerSize)  {  return  Graph (kmerSize);  }

    /** Build a graph from a given bank.
     * \param[in] bank : bank to get the reads from
     * \param[in] fmt : printf-like format for the command line string
     * \return the created graph.
     */
    static Graph  create (bank::IBank* bank, const char* fmt, ...);

    /** Build a graph from user options.
     * \param[in] fmt: printf-like format
     * \return the created graph.
     */
    static Graph  create (const char* fmt, ...);

    /** Build a graph from scratch.
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static Graph  create (tools::misc::IProperties* options)  {  return  Graph (options);  }

    /** Load a graph from some URI.
     * \param[in] uri : the uri to get the graph from
     * \return the loaded graph.
     */
    static Graph  load (const std::string& uri)  {  return  Graph (uri);  }

    /** Get a parser object that knows the user options for building a graph.
     * \return the options parser object.
     */
    static tools::misc::IOptionsParser* getOptionsParser (bool includeMandatory=true, bool enablemphf=false );

    /********************************************************************************/
    template<typename Item, int NB=8>
    class Vector
    {
    public:
        Vector () : _size(0)  {}
        Item& operator[] (size_t idx)  { return (_items[idx]); }
        size_t size()  { return _size; }
        void resize (size_t n)  { _size = n; }

        template<typename Functor>  void iterate (const Functor& f)  { for (size_t i=0; i<_size; i++)  { f(_items[i]); } }

    protected:
        Item   _items[NB];
        size_t _size;
    };

    /********************************************************************************/
    template<typename Item>
    class Iterator : public tools::dp::ISmartIterator<Item>
    {
    public:

        /** */
        Iterator (tools::dp::ISmartIterator<Item>* ref) : _ref(0) { setRef(ref); }

        /** */
        virtual ~Iterator() { setRef(0); }

        /** */
        Iterator (const Iterator<Item>& i) : _ref(0) { setRef(i._ref); }

        /** */
        Iterator& operator= (const Iterator<Item>& i)  {  if (this != &i)   {  setRef (i._ref);  }   return *this;  }

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

    /********************************************************************************/
    /*                               CONSTRUCTORS                                   */
    /********************************************************************************/

    /* Default Constructor.*/
    Graph ();

    /* Copy Constructor.*/
    Graph (const Graph& graph);

    /* Destructor. */
    ~Graph ();

    /** Affectation overload. */
    Graph& operator= (const Graph& graph);

    /**********************************************************************/
    /*                     GLOBAL ITERATOR METHODS                        */
    /**********************************************************************/

    /** Creates an iterator over nodes of the graph.
     * The kind of nodes may depend on the template specialization of this method:
     *      - all nodes for T=Node,
     *      - branching nodes for T=BranchingNode...
     * \return the nodes iterator. */
    template<typename T>
    Graph::Iterator<T> iterator () const;


    /**********************************************************************/
    /*                     ALL NEIGHBORS METHODS                          */
    /**********************************************************************/

    /** Returns a vector of neighbors of the provided node.
     * \param[in] node : the node whose neighbors are wanted
     * \param[in] direction : the direction of the neighbors. If not set, out and in neighbors are computed.
     * \return a vector of the node neighbors (may be empty).
     */
    template <typename T>  Graph::Vector<T> neighbors (const Node& node, Direction direction=DIR_END) const;

    /** Shortcut for 'neighbors' method with direction==DIR_OUTCOMING.
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    template <typename T>  Graph::Vector<T> successors (const Node& node) const;

    /** Shortcut for 'neighbors' method with direction==DIR_INCOMING.
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    template <typename T>  Graph::Vector<T> predecessors (const Node& node) const;

    /** Returns a vector of neighbors of the provided kmer. It has to be understood as the following:
     *  - a node N is built with the kmer, with the strand FORWARD
     *  - a call to 'neighbors<T> (N,          DIR_OUTGOING)' is done; we get v1
     *  - a call to 'neighbors<T> (reverse(N), DIR_OUTGOING)' is done; we get v2
     *  - the result is the concatenation of v1 and v2
     *  \param[in] kmer : the kmer whose neighbors are wanted.
     *  \return a vector of the neighbors (may be empty).
     */
    template <typename T>  Graph::Vector<T> neighbors (const Node::Value& kmer) const;

    /** Returns a set of neighbors for each node iterated with the provided two iterators
     * \param[in] first : beginning of the iteration
     * \param[in] last : end of the iteration
     * \return all the neighbors computed for each iterated node. */
    template <typename T, typename IteratorInput>
    std::set<T> neighbors (IteratorInput first, IteratorInput last) const;

    /** Returns the successors of two nodes, ie with the same transition nucleotide from both nodes.
     * \param[in] node1 : first node
     * \param[in] node2 : sedond node
     * \return the vector of pairs of items as successors
     */
    template <typename T>  Graph::Vector<std::pair<T,T> > successors   (const Node& node1, const Node& node2) const;

    /** Returns the predecessors of two nodes, ie with the same transition nucleotide from both nodes.
     * \param[in] node1 : first node
     * \param[in] node2 : sedond node
     * \return the vector of pairs of items as predecessors
     */
    template <typename T>  Graph::Vector<std::pair<T,T> > predecessors (const Node& node1, const Node& node2) const;

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
    template <typename T>  T neighbor (const Node& source, Direction dir, kmer::Nucleotide nt) const;

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
    template <typename T>  T neighbor (const Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const;

    /** Shortcut for neighbor with dir==DIR_OUTCOMING. */
    template <typename T>  T successor (const Node& source, kmer::Nucleotide nt, bool& exists) const;

    /** Shortcut for neighbor with dir==DIR_OUTCOMING. */
    template <typename T>  T successor (const Node& source, kmer::Nucleotide nt) const;

    /** Shortcut for neighbor with dir==DIR_INCOMING. */
    template <typename T>  T predecessor (const Node& source, kmer::Nucleotide nt, bool& exists) const;

    /** Shortcut for neighbor with dir==DIR_INCOMING. */
    template <typename T>  T predecessor (const Node& source, kmer::Nucleotide nt) const;

    /**********************************************************************/
    /*                      MISC NEIGHBORS METHODS                        */
    /**********************************************************************/

    /** Get the incoming degree of the node.
     * \param[in] node : the node
     * \return the indegree of the node. */
    size_t indegree  (const Node& node) const;

    /** Get the outcoming degree of the node.
     * \param[in] node : the node
     * \return the outdegree of the node. */
    size_t outdegree (const Node& node) const;

    /** Get the degree of the node (either incoming or outcoming).
     * \param[in] node : the node
     * \param[in] dir : direction of the degree
     * \return the degree of the node. */
    size_t degree    (const Node& node, Direction dir) const;

    /** Tells whether or not a transition exists between two given nodes.
     * \param[in] u : first node
     * \param[in] v : second node
     * \return true if such a transition exists, false otherwise. */
    bool isEdge (const Node& u, const Node& v) const;

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
    int simplePathAvance (const Node& node, Direction dir, Edge& output) const;

    /** */
    int simplePathAvance (const Node& node, Direction dir) const;

    /** */
    int simplePathAvance (const Node& node, Direction dir, kmer::Nucleotide& nt) const;

    /** */
    template<typename T> Graph::Iterator<T> simplePath (const Node& node, Direction dir) const;


    /**********************************************************************/
    /*                         NODE METHODS                               */
    /**********************************************************************/

    /** Tells whether or not a node belongs to the graph.
     * \param[in] item : the node
     * \return true if the node belongs to the graph, false otherwise. */
    bool contains (const Node& item) const;

    /** Get the ascii string for the node, according to its strand.
     * \param[in] node: the node to get the string from
     * \return the string representation for the provided node. */
    std::string toString (const Node& node) const;

    /** Tells whether the provided node is branching or not.
     * \param[in] node : the node to be asked
     * \return true if the node is branching, false otherwise. */
    bool isBranching (const Node& node) const;

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
    BranchingNode reverse (const BranchingNode& node) const;

    /** Mutation of a node. */
    Graph::Vector<Node> mutate (const Node& node, size_t idx, int mode=0) const;

    /** Return a nucleotide at position 'idx' of a given node 
     * \param[in] node : the node we want to extract a nucleotide from
     * \param[in] idx : the position of the nucleotide to be extracted
     * \return the wanted nucleotide. */
    kmer::Nucleotide getNT (const Node& node, size_t idx) const;

    /** Return the abundance of a node by querying the perfect hash function 
     * \param[in] node : the node
     * \return the abundance */
    int queryAbundance (const Node& node) const;


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
    std::string toString (const BranchingEdge& edge) const;

    /** Tells whether the provided edge is simple: outdegree(from)==1 and indegree(to)==1
     * \param[in] edge : the edge to be asked
     * \return true if the edge is simple, false otherwise. */
    bool isSimple (const Edge& edge) const;

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
        STATE_BANKCONVERTER_DONE  = (1<<1),
        STATE_SORTING_COUNT_DONE  = (1<<2),
        STATE_BLOOM_DONE          = (1<<3),
        STATE_DEBLOOM_DONE        = (1<<4),
        STATE_BRANCHING_DONE      = (1<<5),
        STATE_MPHF_DONE           = (1<<6)
    };
    typedef int State;
    State getState () const { return _state; }
    bool checkState (StateMask mask) const { return (_state & mask)==mask; }
    State setState   (StateMask mask) { _state |=  mask; return _state; }
    State unsetState (StateMask mask) { _state &= ~mask; return _state; }

private:

    /** Constructor for empty graph.*/
    Graph (size_t kmerSize);

    /** Constructor. Use for Graph creation (ie. DSK + debloom) and filesystem save. */
    Graph (bank::IBank* bank, tools::misc::IProperties* params);

    /** Constructor. Use for Graph creation (ie. DSK + debloom) and filesystem save. */
    Graph (tools::misc::IProperties* params);

    /** Constructor. Use for reading from filesystem. */
    Graph (const std::string& uri);

    /** Kind of storage for the graph. */
    tools::storage::impl::StorageMode_e _storageMode;

    /** Default storage kind. */
    static const tools::storage::impl::StorageMode_e PRODUCT_MODE_DEFAULT = tools::storage::impl::STORAGE_HDF5;

    /** Storage. */
    tools::storage::impl::Storage* _storage;
    void setStorage (tools::storage::impl::Storage* storage)  { SP_SETATTR(storage); }
    tools::storage::impl::Storage& getStorage()                           { return (*_storage); }
    tools::storage::impl::Group&   getGroup  (const std::string name="")  { return getStorage() (name); }

    /** */
    std::string _name;

    /** */
    size_t _kmerSize;

    /** Creation information. */
    tools::misc::impl::Properties _info;

    /** */
    tools::misc::BloomKind       _bloomKind;
    tools::misc::DebloomKind     _debloomKind;
    tools::misc::DebloomImpl     _debloomImpl;
    tools::misc::BranchingKind   _branchingKind;
    tools::misc::MPHFKind        _mphfKind;

    State _state;

    /** Defined as a void* for hiding implementation in cpp file. */
    void* _variant;

    /** */
    Graph::Iterator<Node> getNodes () const;

    /** */
    Graph::Iterator<BranchingNode> getBranchingNodes () const;

    /** */
    Graph::Iterator<Node> getSimpleNodeIterator (const Node& node, Direction dir) const;

    /** */
    Graph::Iterator<Edge> getSimpleEdgeIterator (const Node& node, Direction dir) const;

    /** */
    Graph::Vector<Edge> getEdges (const Node& source, Direction direction) const;

    /** */
    Graph::Vector<std::pair<Node,Node> > getNodesCouple (const Node& node1, const Node& node2, Direction direction) const;

    /** */
    Graph::Vector<std::pair<Edge,Edge> > getEdgesCouple (const Node& node1, const Node& node2, Direction direction) const;

    /** */
    Graph::Vector<Node> getNodes (const Node& source, Direction direction) const;

    /** */
    Graph::Vector<BranchingNode> getBranchingNodeNeighbors (const Node& source, Direction direction) const;

    /** */
    Graph::Vector<BranchingEdge> getBranchingEdgeNeighbors (const Node& source, Direction direction) const;

    /** */
    Graph::Vector<Edge> getEdgeValues (const Node::Value& kmer) const;

    /** */
    Graph::Vector<Node> getNodeValues (const Node::Value& kmer) const;

    /** */
    Graph::Vector<BranchingEdge> getBranchingEdgeValues (const Node::Value& kmer) const;

    /** */
    Graph::Vector<BranchingNode> getBranchingNodeValues (const Node::Value& kmer) const;

    /** */
    Node getNode (const Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const;

    /** Friends. */
    friend struct build_visitor;
};

/********************************************************************************/
/**                           TEMPLATE SPECIALIZATIONS                          */
/********************************************************************************/

template<>   inline Graph::Iterator<Node>          Graph::iterator () const  {  return getNodes ();           }
template<>   inline Graph::Iterator<BranchingNode> Graph::iterator () const  {  return getBranchingNodes ();  }

template <>  inline Graph::Vector<Node> Graph::successors   (const Node& node) const                 {  return getNodes (node, DIR_OUTCOMING); }
template <>  inline Graph::Vector<Node> Graph::predecessors (const Node& node) const                 {  return getNodes (node, DIR_INCOMING);  }
template <>  inline Graph::Vector<Node> Graph::neighbors    (const Node& node, Direction dir) const  {  return getNodes (node, dir);           }
template <>  inline Graph::Vector<Node> Graph::neighbors    (const Node::Value& kmer) const          {  return getNodeValues (kmer);           }

template <>  inline Node Graph::neighbor (const Node& source, Direction dir, kmer::Nucleotide nt) const
{  bool exists=true; return getNode (source, dir, nt, exists);  }

template <>  inline Node Graph::neighbor (const Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{  return getNode (source, dir, nt, exists);  }

template <>  inline Node Graph::successor (const Node& source, kmer::Nucleotide nt) const
{  bool exists=true; return getNode (source, DIR_OUTCOMING, nt, exists);  }

template <>  inline Node Graph::successor (const Node& source, kmer::Nucleotide nt, bool& exists) const
{  return getNode (source, DIR_OUTCOMING, nt, exists);  }

template <>  inline Node Graph::predecessor (const Node& source, kmer::Nucleotide nt) const
{  bool exists=true; return getNode (source, DIR_INCOMING, nt, exists);  }

template <>  inline Node Graph::predecessor (const Node& source, kmer::Nucleotide nt, bool& exists) const
{  return getNode (source, DIR_INCOMING, nt, exists);  }

template <>  inline Graph::Vector<Edge> Graph::successors   (const Node& node) const                 {  return getEdges (node, DIR_OUTCOMING); }
template <>  inline Graph::Vector<Edge> Graph::predecessors (const Node& node) const                 {  return getEdges (node, DIR_INCOMING);  }
template <>  inline Graph::Vector<Edge> Graph::neighbors    (const Node& node, Direction dir) const  {  return getEdges (node, dir);           }
template <>  inline Graph::Vector<Edge> Graph::neighbors    (const Node::Value& kmer) const          {  return getEdgeValues (kmer);           }

/** */
template <>  inline Graph::Vector<std::pair<Edge,Edge> > Graph::successors (const Node& node1, const Node& node2) const
{ return getEdgesCouple (node1, node2, DIR_OUTCOMING); }

template <>  inline Graph::Vector<std::pair<Node,Node> > Graph::successors (const Node& node1, const Node& node2) const
{ return getNodesCouple (node1, node2, DIR_OUTCOMING); }

template <>  inline Graph::Vector<std::pair<Edge,Edge> > Graph::predecessors (const Node& node1, const Node& node2) const
{ return getEdgesCouple (node1, node2, DIR_INCOMING); }

template <>  inline Graph::Vector<std::pair<Node,Node> > Graph::predecessors (const Node& node1, const Node& node2) const
{ return getNodesCouple (node1, node2, DIR_INCOMING); }

/** */
template <>  inline Graph::Vector<BranchingNode> Graph::successors (const Node& node) const
{ return getBranchingNodeNeighbors (node, DIR_OUTCOMING);  }

template <>  inline Graph::Vector<BranchingNode> Graph::predecessors (const Node& node) const
{ return getBranchingNodeNeighbors (node, DIR_INCOMING);  }

/** */
template <>  inline Graph::Vector<BranchingNode> Graph::neighbors (const Node& node, Direction direction) const
{ return getBranchingNodeNeighbors (node, direction);  }

/** */
template <>  inline Graph::Vector<BranchingNode> Graph::neighbors (const Node::Value& kmer) const
{ return getBranchingNodeValues (kmer);  }

/** */
template <> std::set<BranchingNode> Graph::neighbors (std::set<BranchingNode>::iterator first, std::set<BranchingNode>::iterator last) const;

/** */
template <>  inline Graph::Vector<BranchingEdge> Graph::successors (const Node& node) const
{ return getBranchingEdgeNeighbors (node, DIR_OUTCOMING);  }

template <>  inline Graph::Vector<BranchingEdge> Graph::predecessors (const Node& node) const
{ return getBranchingEdgeNeighbors (node, DIR_INCOMING);  }

/** */
template <>  inline Graph::Vector<BranchingEdge> Graph::neighbors (const Node& node, Direction direction) const
{ return getBranchingEdgeNeighbors (node, direction);  }

/** */
template <>  inline Graph::Vector<BranchingEdge> Graph::neighbors (const Node::Value& kmer) const
{ return getBranchingEdgeValues (kmer);  }

/** */
template<> inline Graph::Iterator<Node> Graph::simplePath (const Node& node, Direction dir) const  { return getSimpleNodeIterator(node, dir); }
template<> inline Graph::Iterator<Edge> Graph::simplePath (const Node& node, Direction dir) const  { return getSimpleEdgeIterator(node, dir); }

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
class ProgressGraphIterator : public tools::dp::impl::SubjectIterator<Type>
{
public:
    ProgressGraphIterator (const Graph::Iterator<Type>& items, const char* msg = "compute", size_t divide=100)
        : tools::dp::impl::SubjectIterator<Type> (items.get(), items.size()/divide, new Listener (items.size(), msg)), _size(items.size()) {}

    u_int64_t size () const  { return _size; }

private:
    u_int64_t _size;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
