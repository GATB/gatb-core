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

typedef tools::math::Integer Type;

/** Definition of a Node. */
struct Node
{

    typedef tools::math::Integer Type;

    Node () : abundance(0), strand(kmer::STRAND_FORWARD) {}

    Node (const Type& kmer, kmer::Strand strand=kmer::STRAND_FORWARD, u_int16_t abundance=0) : kmer(kmer), strand(strand), abundance(abundance) {}

    Type         kmer;
    u_int16_t    abundance;
    kmer::Strand strand;

    bool operator== (const Node& other) const  { return kmer == other.kmer; }
    bool operator!= (const Node& other) const  { return kmer != other.kmer; }

    bool operator< (const Node& other) const
    {
        // need to define a strict weak ordering
        if (kmer   != other.kmer)    {  return (kmer   < other.kmer);    }
        if (strand != other.strand)  {  return (strand < other.strand);  }
        return false;
    }

    /** */
    void set (const Type& kmer, const kmer::Strand& strand)
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
        const Type& kmer_from, kmer::Strand strand_from,
        const Type& kmer_to,   kmer::Strand strand_to,
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
    static Graph  create (bank::IBank* bank, tools::misc::IProperties* options=0)
    {
        if (options==0)
        {
            options = new tools::misc::impl::Properties();
            options->add (0, STR_KMER_SIZE, "%d", 27);
            options->add (0, STR_NKS,       "%d", 1);
        }
        return  Graph (bank, options);
    }

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
#if 0
    template<typename Item, int NB=8>
    class Vector
    {
    public:
        Vector () : _size(0)  {}
        Item& operator[] (size_t idx)  { return (_items[idx]); }
        size_t size()  { return _size; }
        void setSize (size_t n)  { _size = n; }

        template<typename Functor>  void iterate (const Functor& f)  { for (size_t i=0; i<_size; i++)  { f(_items[i]); } }

    protected:
        Item   _items[NB];
        size_t _size;
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
    template <typename T>  std::vector<T> successors (const Node& node) const;

    /** */
    template <typename T>  std::vector<T> predecessors (const Node& node) const;

    /** */
    template <typename T>  std::vector<T> neighbors (const Node& node, Direction direction=DIR_END) const;

    /** */
    size_t indegree  (const Node& node) const;
    size_t outdegree (const Node& node) const;
    size_t degree    (const Node& node, Direction dir) const;

    /** */
    bool isEdge (const Node& u, const Node& v) const { return false; }

    /** */
    tools::misc::IProperties& getInfo () const { return (tools::misc::IProperties&)_info; }

    /** */
    std::string toString (const Node& node) const;

    /** */
    void getNearestBranchingRange (const Node& node, Node& begin, Node& end) const;

    /** */
    bool isBranching (const Node& node) const;

    /** */
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
    std::vector<Edge> getEdges (const Node& source, Direction direction) const;

    /** */
    std::vector<Node> getNodes (const Node& source, Direction direction) const;

    /** */
    void executeAlgorithm (tools::misc::impl::Algorithm& algorithm, tools::misc::IProperties* props, tools::misc::IProperties& info);

    /** Friends. */
    template<typename T> friend class GraphFactoryImpl;
};

/********************************************************************************/

template<>   inline Graph::Iterator<Node>          Graph::iterator () const  {  return getNodes ();           }
template<>   inline Graph::Iterator<BranchingNode> Graph::iterator () const  {  return getBranchingNodes ();  }

template <>  inline std::vector<Node> Graph::successors   (const Node& node) const                 {  return getNodes (node, DIR_OUTCOMING); }
template <>  inline std::vector<Node> Graph::predecessors (const Node& node) const                 {  return getNodes (node, DIR_INCOMING);  }
template <>  inline std::vector<Node> Graph::neighbors    (const Node& node, Direction dir) const  {  return getNodes (node, dir);           }

/********************************************************************************/

template <>  inline std::vector<Edge> Graph::successors   (const Node& node) const                 {  return getEdges (node, DIR_OUTCOMING); }
template <>  inline std::vector<Edge> Graph::predecessors (const Node& node) const                 {  return getEdges (node, DIR_INCOMING);  }
template <>  inline std::vector<Edge> Graph::neighbors    (const Node& node, Direction dir) const  {  return getEdges (node, dir);           }

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
