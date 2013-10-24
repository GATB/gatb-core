/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file IGraph.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface for DeBruijn graph
 */

#ifndef _GATB_CORE_DEBRUIJN_IGRAPH_HPP_
#define _GATB_CORE_DEBRUIJN_IGRAPH_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/kmer/api/IModel.hpp>
#include <sstream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
/********************************************************************************/

enum Direction
{
    DIR_OUTCOMING = 0,
    DIR_INCOMING  = 1
};

/********************************************************************************/

typedef tools::math::Integer Type;

/** Forward declaration. */
class IGraph;
template<typename Item, int NB> class GraphItemSet;

/********************************************************************************/

/** Definition of a graph item: a node or a edge (so far). */
class GraphItem
{
public:

    const IGraph* getGraph() const  { return graph; }

protected:

    /** */
    const IGraph* graph;
    void setGraph (const IGraph* g)  { this->graph = g; }

    template<typename Item, int NB>  friend class GraphItemSet;
};

/********************************************************************************/

template<typename Item, int NB=8>
class GraphItemSet
{
public:
    GraphItemSet (const IGraph* graph) : _size(0)
    {
        for (size_t i=0; i<NB; i++)  { _items[i].setGraph (graph); }
    }

    ~GraphItemSet () {}

    Item& operator[] (size_t idx)  { return (_items[idx]); }

    size_t size()  { return _size; }

    void setSize (size_t n)  { _size = n; }

protected:
    Item   _items[NB];
    size_t _size;
};

/********************************************************************************/

/** Definition of a Node. */
class Node : public GraphItem
{
public:
    kmer::Kmer<Type>  kmer;
    kmer::Strand      strand;

//private:

    /** */
    void set (const IGraph* graph, const Type& k, const kmer::Strand& s) { this->graph=graph; kmer.value=k; strand=s; }

    friend class Edge;
    friend class Graph;
};

/********************************************************************************/

/** Definition of an Edge. */
class Edge : public GraphItem
{
public:

    Node       from;
    Node       to;
    kmer::Nucleotide nt;
    Direction  direction;

//private:

    void set (
        const IGraph* graph,
        const Type& kmer_from, kmer::Strand strand_from,
        const Type& kmer_to,   kmer::Strand strand_to,
        kmer::Nucleotide n, Direction dir
    )
    {
        from.set (graph, kmer_from, strand_from);
        to.set   (graph, kmer_to,   strand_to);
        nt = n;
        direction = dir;
    }

    friend class Graph;
};

/********************************************************************************/
class NodeSet : public GraphItemSet<Node>
{
public: NodeSet (const IGraph& graph) : GraphItemSet<Node>(&graph) {}
};

class EdgeSet : public GraphItemSet<Edge>
{
public: EdgeSet (const IGraph& graph) : GraphItemSet<Edge>(&graph) {}
};

/********************************************************************************/
class INodeIterator : public tools::dp::Iterator<Node>
{
public:

    /** */
    virtual ~INodeIterator() {}

    /** */
    virtual u_int64_t getNbItems () const = 0;
};

/********************************************************************************/

template<class Listener>
class ProgressIterator : public tools::dp::impl::SubjectIterator<Node>
{
public:
    ProgressIterator (INodeIterator* nodes, const char* msg = "compute")
        : tools::dp::impl::SubjectIterator<Node> (nodes, nodes->getNbItems()/100, new Listener (nodes->getNbItems(), msg)) {}
};

/********************************************************************************/

/**
 * Inspired from http://jung.sourceforge.net/doc/api/index.html
 */
class IGraph : public tools::dp::SmartPointer
{
public:

    /** Destructor. */
    virtual ~IGraph() {}

    /** From Container interface. */
    virtual bool contains (const Node& item) const = 0;

    /** Creates an iterator over all nodes of the graph.
     * \return the all nodes iterator. */
    virtual INodeIterator* nodes () const = 0;

    /** */
    virtual size_t getOutEdges (const Node& node, EdgeSet& edges) const =  0;

    /** */
    virtual size_t getInEdges (const Node& node, EdgeSet& edges) const =  0;

    /** */
    virtual size_t getSuccessors (const Node& node, NodeSet& nodes) const = 0;

    /** */
    virtual size_t getPredecessors (const Node& node, NodeSet& nodes) const = 0;

    /** */
    virtual bool isEdge (const Node& u, const Node& v) const = 0;

    /** */
    virtual tools::misc::IProperties& getInfo () const = 0;

    /** */
    virtual std::string toString (const Node& node, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const = 0;

    /** Remove physically the graph. */
    virtual void remove () = 0;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IGRAPH_HPP_ */
