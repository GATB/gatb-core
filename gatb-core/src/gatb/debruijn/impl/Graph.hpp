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
    DIR_OUTCOMING = 0,
    DIR_INCOMING  = 1
};

/********************************************************************************/

typedef tools::math::Integer Type;

/** Definition of a Node. */
struct Node
{
    Type         kmer;
    u_int16_t    abundance;
    kmer::Strand strand;

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
};

/********************************************************************************/

class Graph
{
public:

    /********************************************************************************/
    /** Build a graph from a given bank.
     * \param[in] bank : bank to get the reads from
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static Graph  create (bank::IBank* bank, tools::misc::IProperties* options)  {  return  Graph (bank, options);  }

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
    template<typename Item, int NB=8>
    class Vector
    {
    public:
        Vector () : _size(0)  {}
        Item& operator[] (size_t idx)  { return (_items[idx]); }
        size_t size()  { return _size; }
        void setSize (size_t n)  { _size = n; }
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

        template<class Listener> friend class ProgressIterator;
    };

    /** Destructor. */
    ~Graph ();

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

    /** */
    bool isEdge (const Node& u, const Node& v) const { return false; }

    /** */
    tools::misc::IProperties& getInfo () const { return (tools::misc::IProperties&)_info; }

    /** */
    std::string toString (const Node& node, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const;

    /** */
    void getNearestBranchingRange (const Node& node, Node& begin, Node& end) const;

    /** */
    bool isBranching (const Node& node) const;

private:

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
    void executeAlgorithm (tools::misc::impl::Algorithm& algorithm, tools::misc::IProperties* props, tools::misc::IProperties& info);

    /** Friends. */
    template<typename T> friend class GraphFactoryImpl;
};

/********************************************************************************/

template<>   inline Graph::Iterator<Node>          Graph::iterator () const  {  return getNodes ();           }
template<>   inline Graph::Iterator<BranchingNode> Graph::iterator () const  {  return getBranchingNodes ();  }

template <>  inline Graph::Vector<Node> Graph::successors   (const Node& node) const  {  return getNodes (node, DIR_OUTCOMING); }
template <>  inline Graph::Vector<Node> Graph::predecessors (const Node& node) const  {  return getNodes (node, DIR_INCOMING);  }

/********************************************************************************/

template <>  inline Graph::Vector<Edge> Graph::successors   (const Node& node) const  {  return getEdges (node, DIR_OUTCOMING); }
template <>  inline Graph::Vector<Edge> Graph::predecessors (const Node& node) const  {  return getEdges (node, DIR_INCOMING);  }

/********************************************************************************/

template<class Listener>
class ProgressIterator : public tools::dp::impl::SubjectIterator<Node>
{
public:
    ProgressIterator (const Graph::Iterator<Node>& nodes, const char* msg = "compute")
        : tools::dp::impl::SubjectIterator<Node> (nodes.get(), nodes.getNbItems()/100, new Listener (nodes.getNbItems(), msg)) {}
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
