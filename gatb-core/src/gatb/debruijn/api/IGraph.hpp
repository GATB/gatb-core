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

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <sstream>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for DeBruijn graph management. */
namespace debruijn  {
/********************************************************************************/

enum Strand
{
    STRAND_FORWARD = (1<<0),
    STRAND_REVCOMP = (1<<1),
    STRAND_ALL     = STRAND_FORWARD + STRAND_REVCOMP
};

enum Nucleotide
{
    NUCL_A   = (1<<0),
    NUCL_C   = (1<<1),
    NUCL_T   = (1<<2),
    NUCL_G   = (1<<3),
    NUCL_ALL = NUCL_A + NUCL_C + NUCL_T + NUCL_G
};

enum Direction
{
    DIR_OUTCOMING = 0,
    DIR_INCOMING  = 1
};

/********************************************************************************/

/** Forward declaration. */
template<typename T> class IGraph;
template<typename T> class Graph;
template<typename T,typename Item> class GraphItemSet;

/********************************************************************************/

/** Definition of a graph item: a node or a edge (so far). */
template<typename T>
class GraphItem
{
public:

    /** Constructor. */
    GraphItem (IGraph<T>* graph) : _graph(graph)  {}

protected:

    /** */
    IGraph<T>* _graph;
};

/********************************************************************************/

template<typename T, typename Item>
class GraphItemSet
{
public:
    GraphItemSet () : _size(0)  { }

    ~GraphItemSet ()  {  for (size_t i=0; i<8; i++)  {  delete _items[i]; }  }


    Item& operator[] (size_t idx)  { return *(_items[idx]); }

    size_t size()  { return _size; }

    void setSize (size_t n)  { _size = n; }

protected:
    Item*  _items[8];
    size_t _size;
};

/********************************************************************************/


/** Definition of a Node. */
template<typename T>
class Node : public GraphItem<T>
{
public:
    T         kmer;
    Strand    strand;

    Node (IGraph<T>* graph=0) : GraphItem<T> (graph), kmer(0), strand(STRAND_FORWARD)  {}

    /** */
    void set (const T& k, const Strand& s) { kmer=k; strand=s; }

    /** */
    Node reverse ()
    {
        Node res;
        res.kmer = this->_graph->getModel().revcomp(kmer);
        res.strand = strand == STRAND_FORWARD ? STRAND_REVCOMP : STRAND_FORWARD;
        return res;
    }

    /** For debug mainly. */
    std::string toString (Strand strand = STRAND_ALL)
    {
        std::stringstream ss;

        if (strand == STRAND_ALL || this->strand == strand)
        {
            ss << "[ " << this->_graph->getModel().toString (kmer) <<  "  strand=" << (this->strand==STRAND_FORWARD ? "FORWARD" : "REVCOMP") << "]";
        }
        else
        {
            T reverse = this->_graph->getModel().reverse (kmer);
            ss << "[ " << this->_graph->getModel().toString (reverse) <<  "  strand=" << (this->strand==STRAND_FORWARD ? "REVCOMP" : "FORWARD") << "]";
        }
        return ss.str();
    }

    /** Returns a new node with the same kmer but with the opposite strand. */
    friend Node operator-(const Node<T>& in)
    {
        Node res;
        res.kmer = in.kmer;
        res.strand = in.strand == STRAND_FORWARD ? STRAND_REVCOMP : STRAND_FORWARD;
        return res;
    }

    friend class GraphItemSet<T,Node<T> >;
};

/********************************************************************************/

/** Definition of an Edge. */
template<typename T>
class Edge : public GraphItem<T>
{
public:

    Edge (IGraph<T>* graph) : GraphItem<T> (graph), from(graph), to(graph), nt(NUCL_ALL), direction(DIR_OUTCOMING)  {}

    Node<T>    from;
    Node<T>    to;
    Nucleotide nt;
    Direction  direction;

    /** */
    void setGraph (IGraph<T>* graph)
    {
        this->_graph = graph;
        from.setGraph (graph);
        to.setGraph (graph);
    }

    /** */
    void set (
        const T& kmer_from, Strand strand_from,
        const T& kmer_to,   Strand strand_to,
        Nucleotide n, Direction dir
    )
    {
        from.set (kmer_from, strand_from);
        to.set   (kmer_to,   strand_to);
        nt = n;
        direction = dir;
    }

    /** */
    Edge& operator= (const Edge& e)
    {
        if (this != &e)   { from = e.from;  to = e.to;  nt = e.nt;  direction = e.direction; }
        return *this;
    }

    std::string toString (Strand strand = STRAND_ALL)
    {
        std::stringstream ss;
        static char nucl[] = {'A', 'C', 'T', 'G'};
        ss << "[EDGE: "
           << "nt=" << nucl[nt] << "  "
           << "dir=" << (direction==DIR_INCOMING ? "in" : "out") << "  "
           << "from " << from.toString(strand)   << " "
           << "to "   << to.toString(strand)     << " "
           << "]";
        return ss.str();
    }

    friend class GraphItemSet<T,Edge<T> >;
};


/********************************************************************************/
template<typename T>  class NodeSet : public GraphItemSet<T, Node<T> >
{
public:
    NodeSet (Graph<T>& graph)
    {
        /** For each item, we set the reference on the graph. */
        for (size_t i=0; i<8; i++)  {  this->_items[i] = new Node<T> (graph._ref);  }
     }
};

template<typename T>  class EdgeSet : public GraphItemSet<T, Edge<T> >
{
public:
    EdgeSet (Graph<T>& graph)
    {
        /** For each item, we set the reference on the graph. */
        for (size_t i=0; i<8; i++)  {  this->_items[i] = new Edge<T> (graph._ref);  }
    }
};

/********************************************************************************/
template<typename T>
class INodeIterator : public tools::dp::Iterator<Node<T> >
{
};

/********************************************************************************/
template<typename T>
class NodeIterator : public INodeIterator<T>
{
public:

    /** Constructor. */
    NodeIterator (INodeIterator<T>* ref) : _ref(0)  { setRef(ref); }

    /** Destructor. */
    ~NodeIterator ()  { setRef(0); }

    /** \copydoc  Iterator::first */
    void first() { _ref->first(); }

    /** \copydoc  Iterator::next */
    void next()  { _ref->next(); }

    /** \copydoc  Iterator::isDone */
    bool isDone() { return _ref->isDone();  }

    /** \copydoc  Iterator::item */
    Node<T>& item ()  {  return _ref->item();  }

    /** \copydoc  Iterator::setItem */
    void setItem (Node<T>& i)  { _ref->setItem(i); }

protected:

    INodeIterator<T>* _ref;
    void setRef (INodeIterator<T>* ref)  { SP_SETATTR(ref); }
};

/********************************************************************************/

/**
 * Inspired from http://jung.sourceforge.net/doc/api/index.html
 */
template<typename T>
class IGraph : public tools::dp::SmartPointer
{
public:

    /** Destructor. */
    virtual ~IGraph() {}

    /** From Container interface. */
    virtual bool contains (const Node<T>& item) = 0;

    /** Creates an iterator over all nodes of the graph.
     * \return the all nodes iterator. */
    virtual INodeIterator<T>* nodes () = 0;

    /** */
    virtual size_t getOutEdges (const Node<T>& node, EdgeSet<T>& edges) =  0;

    /** */
    virtual size_t getInEdges (const Node<T>& node, EdgeSet<T>& edges) =  0;

    /** */
    virtual size_t getSuccessors (const Node<T>& node, NodeSet<T>& nodes) = 0;

    /** */
    virtual size_t getPredecessors (const Node<T>& node, NodeSet<T>& nodes) = 0;

    /** */
    virtual bool isEdge (const Node<T>& u, const Node<T>& v) = 0;

    /** */
    virtual kmer::impl::Model<T>& getModel ()  = 0;
};

/********************************************************************************/

template<typename T>
class Graph
{
public:

    /** Constructor. */
    Graph (IGraph<T>* ref) : _ref(0)  { setRef(ref); }

    /** Destructor. */
    ~Graph () { setRef(0); }

    /** From Container interface. */
    bool contains (const Node<T>& item)  { return _ref->contains (item); }

    /** Creates an iterator over all nodes of the graph.
     * \return the all nodes iterator. */
    NodeIterator<T> nodes () { return NodeIterator<T> (_ref->nodes ()); }

    /** */
    size_t getOutEdges (const Node<T>& node, EdgeSet<T>& edges)     { return _ref->getOutEdges (node, edges); }

    /** */
    size_t getInEdges (const Node<T>& node, EdgeSet<T>& edges)      { return _ref->getInEdges (node, edges); }

    /** */
    size_t getSuccessors (const Node<T>& node, NodeSet<T>& nodes)   { return _ref->getSuccessors (node, nodes); }

    /** */
    size_t getPredecessors (const Node<T>& node, NodeSet<T>& nodes) { return _ref->getPredecessors (node, nodes); }

    /** */
    kmer::impl::Model<T>& getModel ()  { return _ref->getModel(); }

//private:

    IGraph<T>* _ref;
    void setRef (IGraph<T>* ref)  { SP_SETATTR(ref); }
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IGRAPH_HPP_ */
