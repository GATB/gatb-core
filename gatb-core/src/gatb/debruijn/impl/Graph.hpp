/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file GraphBasic.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HPP_

/********************************************************************************/

//#include <gatb/debruijn/api/IGraph.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/kmer/api/IModel.hpp>

#include <gatb/bank/api/IBank.hpp>

#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <vector>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {

/********************************************************************************/

#if 0
#define ProductFactoryLocal tools::collections::impl::ProductFileFactory
#else
#define ProductFactoryLocal tools::collections::impl::ProductHDF5Factory
#endif

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

/********************************************************************************/

//typedef tools::math::Integer Type;

/** Definition of a Node. */
struct Node
{

    typedef tools::math::Integer Value;

    Node () : abundance(0), strand(kmer::STRAND_FORWARD) {}

    Node (const Node::Value& kmer, kmer::Strand strand=kmer::STRAND_FORWARD, u_int16_t abundance=0) : kmer(kmer), strand(strand), abundance(abundance) {}

    Node::Value  kmer;
    u_int16_t    abundance;
    kmer::Strand strand;

    bool operator== (const Node& other) const  { return kmer == other.kmer; }
    bool operator!= (const Node& other) const  { return kmer != other.kmer; }

    bool operator< (const Node& other) const  { return (kmer   < other.kmer); }

    /** */
    void set (const Node::Value& kmer, const kmer::Strand& strand)
    {
        this->kmer      = kmer;
        this->strand    = strand;
    }
};

/********************************************************************************/

struct BranchingNode : Node
{
};

/********************************************************************************/

/** Definition of an Edge. */
struct Edge
{
    Node             from;
    Node             to;
    kmer::Nucleotide nt;
    Direction        direction;

    /** */
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

#if 0
    /** */
    Edge reverse() const
    {
        Edge result;
        result.set (
            to.kmer,   to.strand,
            from.kmer, from.strand,
            nt,
            direction==DIR_OUTCOMING ? DIR_INCOMING : DIR_OUTCOMING
        );
        return result;
    }
#endif
};

/********************************************************************************/

class Graph
{
public:

    /********************************************************************************/
    /** Build an empty graph.
     * \param[in] kmerSize: kmer size
     * \return the created graph.
     */
    static Graph  create (size_t kmerSize=27)  {  return  Graph (kmerSize);  }

    /** Build a graph from a given bank.
     * \param[in] bank : bank to get the reads from
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static Graph  create (bank::IBank* bank, const char* fmt, ...);

    /** Build a graph from scratch.
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static Graph  create (tools::misc::IProperties* options)  {  return  Graph (options);  }

    /** Load a graph from some URI.
     * \parm[in] uri : the uri to get the graph from
     * \return the loaded graph.
     */
    static Graph  load (const std::string& uri)  {  return  Graph (uri);  }

    /********************************************************************************/
#if 1
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
#else
    /** Not optimal since std::vector uses dynamic allocation, although we know we can't have more that 8 items.*/
    template<typename Item>
    class Vector : public std::vector<Item>
    {
    public:
        Vector () : std::vector<Item>(8) {}
        template<typename Functor>  void iterate (const Functor& f)  { for (size_t i=0; i<this->size(); i++)  { f((*this)[i]); } }
    };
#endif
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
        u_int64_t getNbItems () const  { return _ref->getNbItems(); }

        /** */
        tools::dp::ISmartIterator<Item>* get()  const { return _ref; }

    private:

        tools::dp::ISmartIterator<Item>* _ref;
        void setRef (tools::dp::ISmartIterator<Item>* ref)  { SP_SETATTR(ref); }
    };

    /********************************************************************************/

    /** Default Constructor.*/
    Graph ();

    /** Copy Constructor.*/
    Graph (const Graph& graph);

    /** Destructor. */
    ~Graph ();

    /** */
    Graph& operator= (const Graph& graph);

    /** */
    void remove ();

    /** From Container interface. */
    bool contains (const Node& item) const;

    /** Creates an iterator over all nodes of the graph.
     * \return the all nodes iterator. */
    template<typename T>
    Graph::Iterator<T> iterator () const;

    /** */
    template <typename T>  Graph::Vector<T> successors (const Node& node) const;

    /** */
    template <typename T>  Graph::Vector<T> predecessors (const Node& node) const;

    /** Returns a vector of neighbors of the provided node.
     * \param[in] node : the node whose neighbors are wanted
     * \param[in] direction : the direction of the neighbors.
     * \return a vector of the node neighbors (may be empty). */
    template <typename T>  Graph::Vector<T> neighbors (const Node& node, Direction direction=DIR_END) const;

    /** */
    template <typename T>  Graph::Vector<T> neighbors (const Node::Value& kmer) const;

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
     * \param[in] direction : direction of the degree
     * \return the degree of the node. */
    size_t degree    (const Node& node, Direction dir) const;

    /** */
    bool isEdge (const Node& u, const Node& v) const { return false; }

    /** Get information about the graph (gathered during its creation).
     * \return a property object holding graph information. */
    tools::misc::IProperties& getInfo () const { return (tools::misc::IProperties&)_info; }

    /** Get the ascii string for the node, according to its strand.
     * \param[in] node: the node to get the string from
     * \return the string representation for the provided node. */
    std::string toString (const Node& node) const;

    /** */
    void getNearestBranchingRange (const Node& node, Node& begin, Node& end) const;

    /** Tells whether the provided node is branching or not.
     * \param[in] node : the node to be asked
     * \return true if the node is branching, false otherwise. */
    bool isBranching (const Node& node) const;

    /** Get the size of the kmers.
     * \return the kmer size. */
    size_t getKmerSize() const { return _kmerSize; }

    /** */
    Node getNode (const tools::misc::Data& data, size_t offset=0) const;

    /** */
    Node reverse (const Node& node) const;

    /** */
    Edge reverse (const Edge& edge) const;

    /**********************************************************************/
    /*                         DEBUG METHODS                              */
    /**********************************************************************/
    /** */
    std::string debugString (const Node& node, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const;

    /** */
    std::string debugString (const Edge& edge, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const;

private:

    /** Constructor for empty graph.*/
    Graph (size_t kmerSize);

    /** Constructor. Use for Graph creation (ie. DSK + debloom) and filesystem save. */
    Graph (bank::IBank* bank, tools::misc::IProperties* params);

    /** Constructor. Use for Graph creation (ie. DSK + debloom) and filesystem save. */
    Graph (tools::misc::IProperties* params);

    /** Constructor. Use for reading from filesystem. */
    Graph (const std::string& uri);

    /** Product. */
    tools::collections::impl::Product<ProductFactoryLocal>* _product;
    void setProduct (tools::collections::impl::Product<ProductFactoryLocal>* product)  { SP_SETATTR(product); }
    tools::collections::impl::Group<ProductFactoryLocal>& getProduct(const std::string name="")  { return (*_product) (name); }

    /** */
    size_t _kmerSize;

    /** Creation information. */
    tools::misc::impl::Properties _info;

    /** Defined as a void* for hiding implementation in cpp file. */
    void* _variant;

    /** */
    Graph::Iterator<Node> getNodes () const;

    /** */
    Graph::Iterator<BranchingNode> getBranchingNodes () const;

    /** */
    Graph::Vector<Edge> getEdges (const Node& source, Direction direction) const;

    /** */
    Graph::Vector<Node> getNodes (const Node& source, Direction direction) const;

    /** */
    Graph::Vector<Edge> getEdgeValues (const Node::Value& kmer) const;

    /** */
    Graph::Vector<Node> getNodeValues (const Node::Value& kmer) const;

    /** */
    void executeAlgorithm (tools::misc::impl::Algorithm& algorithm, tools::misc::IProperties* props, tools::misc::IProperties& info);

    /** Friends. */
    template<typename T> friend class GraphFactoryImpl;
};

/********************************************************************************/

template<>   inline Graph::Iterator<Node>          Graph::iterator () const  {  return getNodes ();           }
template<>   inline Graph::Iterator<BranchingNode> Graph::iterator () const  {  return getBranchingNodes ();  }

template <>  inline Graph::Vector<Node> Graph::successors   (const Node& node) const                 {  return getNodes (node, DIR_OUTCOMING); }
template <>  inline Graph::Vector<Node> Graph::predecessors (const Node& node) const                 {  return getNodes (node, DIR_INCOMING);  }
template <>  inline Graph::Vector<Node> Graph::neighbors    (const Node& node, Direction dir) const  {  return getNodes (node, dir);           }
template <>  inline Graph::Vector<Node> Graph::neighbors    (const Node::Value& kmer) const          {  return getNodeValues (kmer);           }

/********************************************************************************/

template <>  inline Graph::Vector<Edge> Graph::successors   (const Node& node) const                 {  return getEdges (node, DIR_OUTCOMING); }
template <>  inline Graph::Vector<Edge> Graph::predecessors (const Node& node) const                 {  return getEdges (node, DIR_INCOMING);  }
template <>  inline Graph::Vector<Edge> Graph::neighbors    (const Node& node, Direction dir) const  {  return getEdges (node, dir);           }
template <>  inline Graph::Vector<Edge> Graph::neighbors    (const Node::Value& kmer) const          {  return getEdgeValues (kmer);           }

/********************************************************************************/

template<class Type, class Listener>
class ProgressIterator : public tools::dp::impl::SubjectIterator<Type>
{
public:
    ProgressIterator (const Graph::Iterator<Type>& items, const char* msg = "compute")
        : tools::dp::impl::SubjectIterator<Type> (items.get(), items.getNbItems()/100, new Listener (items.getNbItems(), msg)) {}
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
