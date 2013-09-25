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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_

/********************************************************************************/

#include <gatb/debruijn/api/IGraph.hpp>
#include <gatb/tools/designpattern/api/SmartPointer.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/collections/impl/ContainerSet.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Algorithm.hpp>

#include <gatb/tools/collections/impl/Product.hpp>
#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

#include <math.h>

#if 0
#define ProductFactoryLocal tools::collections::impl::ProductFileFactory
#else
#define ProductFactoryLocal tools::collections::impl::ProductHDF5Factory
#endif

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {

/********************************************************************************/

template<typename T>
class GraphBasic : public IGraph<T>
{
public:

    /** Constructor. */
    GraphBasic (
        tools::collections::Iterable<kmer::Kmer<T> >* solidKmers,
        tools::collections::Iterable<T>*              cFPKmers,
        size_t kmerSize
    );

    /** Constructor. */
    GraphBasic (tools::misc::IProperties* options);

    /** Constructor. */
    GraphBasic (const std::string& uri);

    /** Destructor. */
    ~GraphBasic ();

    /** */
    void remove () { _product->remove (); }

    /** */
    void build (bank::IBank* bank);

    /** */
    bool contains (const Node<T>& item)  { return contains (item.kmer.value); }

    /** */
    INodeIterator<T>* nodes ()  { return new NodeIteratorBasic (this, _solidKmersIterable->iterator(), _solidKmersIterable->getNbItems()); }

    /** */
    bool isEdge (const Node<T>& u, const Node<T>& v)  { return false; }

    /** */
    size_t getOutEdges     (const Node<T>& node, EdgeSet<T>& edges) {  return getEdges (node, edges, DIR_OUTCOMING); }

    /** */
    size_t getInEdges      (const Node<T>& node, EdgeSet<T>& edges) {  return getEdges (node, edges, DIR_INCOMING); }

    /** */
    size_t getSuccessors   (const Node<T>& node, NodeSet<T>& nodes) { return getNodes  (node, nodes, DIR_OUTCOMING); }

    /** */
    size_t getPredecessors (const Node<T>& node, NodeSet<T>& nodes) { return getNodes  (node, nodes, DIR_INCOMING); }

    kmer::impl::Model<T>& getModel ()  { return *_model; }

    tools::misc::IProperties* getInfo ()  { return &_info; }

private:

    /** We define states for the graph. */
    enum State_e
    {
        INIT,
        READY
    };
    State_e _state;
    State_e getState () const { return _state; }
    void setState (State_e state)  { _state = state; }

    /** Shortcuts. */
    size_t _kmerSize;
    size_t _nks;
    T      _mask;
    size_t _span;

    tools::collections::impl::Product<ProductFactoryLocal>* _product;
    void setProduct (tools::collections::impl::Product<ProductFactoryLocal>* product)  { SP_SETATTR(product); }
    tools::collections::impl::Group<ProductFactoryLocal>& getProduct(const std::string name="")  { return (*_product) (name); }

    tools::collections::Iterable<kmer::Kmer<T> >* _solidKmersIterable;
    void setSolidKmersIterable (tools::collections::Iterable<kmer::Kmer<T> >* solidKmersIterable)  { SP_SETATTR(solidKmersIterable); }

    tools::collections::Container<T>* _cFPset;
    void setcFPset (tools::collections::Container<T>* cFPset)  { SP_SETATTR(cFPset); }

    tools::collections::impl::Bloom<T>* _bloom;
    void setBloom (tools::collections::impl::Bloom<T>* bloom)  { SP_SETATTR(bloom); }

    tools::misc::IProperties* _props;
    void setProps (tools::misc::IProperties* props)  { SP_SETATTR(props); }

    kmer::impl::Model<T>* _model;
    void setModel (kmer::impl::Model<T>* model)
    {
        if (_model != 0)  { delete _model; }
        _model = model;
    }

    tools::misc::impl::Properties _info;

    /** */
    void configure (
        size_t kmerSize,
        tools::collections::Iterable<kmer::Kmer<T> >* solidIterable,
        tools::collections::Iterable<T>*              cFPKmers,
        tools::collections::Iterable<tools::math::NativeInt8>* bloomArray=0
    );

    /** */
    void executeAlgorithm (gatb::core::tools::misc::impl::Algorithm& algorithm);

    /** */
    bool contains (const T& item)  { return _bloom->contains (item) && !_cFPset->contains(item); }

    /** */
    size_t getEdges (const Node<T>& source, EdgeSet<T>& edges, Direction direction);

    /** */
    size_t getNodes (const Node<T>& source, NodeSet<T>& nodes, Direction direction);

    /**********************************************************************/
    class NodeIteratorBasic : public INodeIterator<T>
    {
    public:

        NodeIteratorBasic (GraphBasic* graph, tools::dp::Iterator<kmer::Kmer<T> >* ref, u_int64_t nbItems)
            : _ref(0), _graph(graph), _isDone(true), _nbItems(nbItems)
        {
            setRef(ref);
            (this->_item)->setGraph (_graph);
            _ref->setItem (this->_item->kmer);
        }

        ~NodeIteratorBasic ()  { setRef(0);   }

        /** \copydoc  Iterator::first */
        void first()
        {
            _ref->first();
            _isDone = _ref->isDone();
        }

        /** \copydoc  Iterator::next */
        void next()
        {
            _ref->next();
            _isDone = _ref->isDone();
        }

        /** \copydoc  Iterator::isDone */
        bool isDone() { return _isDone;  }

        /** \copydoc  Iterator::item */
        Node<T>& item ()  {  return *(this->_item);  }

        /** */
        void setItem (Node<T>& i)
        {
            /** We set the node item to be set for the current iterator. */
            this->_item = &i;

            /** We set the graph for the current node. */
            // (this->_item)->setGraph (_graph);

            /** We set the kmer item to be set for the kmer iterator. */
            _ref->setItem (this->_item->kmer);
        }

        /** */
        u_int64_t getNbItems () const { return _nbItems; }

    private:
        tools::dp::Iterator<kmer::Kmer<T> >* _ref;
        void setRef (tools::dp::Iterator<kmer::Kmer<T> >* ref)  { SP_SETATTR(ref); }

        GraphBasic* _graph;

        bool _isDone;

        u_int64_t _nbItems;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
