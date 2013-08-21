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
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/collections/impl/ContainerSet.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/system/impl/System.hpp>

#include <math.h>

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
    GraphBasic (tools::misc::IProperties* props)
        : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0)
    {
        setProps (props);

        setNodesIterable (new tools::collections::impl::IterableFile<T> (props->getStr (STR_KMER_SOLID)));

        setModel (new kmer::impl::Model<T> (props->getInt(STR_KMER_SIZE)));

        /** We compute the bloom size. */
        double lg2 = log(2);
        float NBITS_PER_KMER = log (16*props->getInt(STR_KMER_SIZE)*(lg2*lg2))/(lg2*lg2);
        u_int64_t bloomSize =  (u_int64_t) (system::impl::System::file().getSize (props->getStr (STR_KMER_SOLID)) / sizeof (T) * NBITS_PER_KMER);

        setBloom (kmer::impl::BloomBuilder<T> (_nodesIterable->iterator(), bloomSize, 7).build ());

        if (props->get(STR_KMER_CFP) != 0)
        {
            setcFPset (new tools::collections::impl::ContainerSet<T> (
                new tools::collections::impl::IteratorFile<T> (props->getStr(STR_KMER_CFP))
            ));
        }

        _mask = _model->getMask ();
        _span = _model->getSpan ();
    }

    /** Constructor. */
    GraphBasic (
        tools::collections::Iterable<T>* solidKmers,
        tools::collections::Iterable<T>* cFPKmers,
        size_t kmerSize
    )
        : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0)
    {
        setNodesIterable (solidKmers);

        setModel (new kmer::impl::Model<T> (kmerSize));

        /** We compute the bloom size. */
        double lg2 = log(2);
        float NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
        u_int64_t bloomSize =  (u_int64_t) (solidKmers->getNbItems() * NBITS_PER_KMER);

        setBloom (kmer::impl::BloomBuilder<T> (_nodesIterable->iterator(), bloomSize, 7).build ());

        setcFPset (new tools::collections::impl::ContainerSet<T> (cFPKmers->iterator()));

        _mask = _model->getMask ();
        _span = _model->getSpan ();
    }

    /** Destructor. */
    ~GraphBasic ()
    {
        setProps         (0);
        setNodesIterable (0);
        setcFPset        (0);
        setModel         (0);
        setBloom         (0);
    }

    /** */
    bool contains (const Node<T>& item)  { return contains (item.kmer); }

    /** */
    INodeIterator<T>* nodes ()  { return new NodeIteratorBasic (this, _nodesIterable->iterator()); }

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

private:

    T      _mask;
    size_t _span;


    tools::collections::Iterable<T>* _nodesIterable;
    void setNodesIterable (tools::collections::Iterable<T>* nodesIterable)
    {
        //if (_nodesIterable != 0)  { delete _nodesIterable; }
        _nodesIterable = nodesIterable;
    }

    tools::collections::Container<T>* _cFPset;
    void setcFPset (tools::collections::Container<T>* cFPset)
    {
        //if (_cFPset != 0)  { delete _cFPset; }
        _cFPset = cFPset;
    }

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

    /** */
    bool contains (const T& item)  { return _bloom->contains (item) && !_cFPset->contains(item); }


    /**********************************************************************/
    size_t getEdges (const Node<T>& source, EdgeSet<T>& edges, Direction direction)
    {
        size_t idx = 0;

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        T graine = ((source.strand == STRAND_FORWARD) ?  source.kmer : revcomp (source.kmer, _span));

        if (direction == DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ( (graine << 2 )  + nt) & _mask;
                T reverse = revcomp (forward, _span);

                if (forward < reverse)
                {
                    if (this->contains (forward))
                    {
                        edges[idx++].set (source.kmer, source.strand, forward, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
                else
                {
                    if (this->contains (reverse))
                    {
                        edges[idx++].set (source.kmer, source.strand, reverse, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
            }
        }
        else if (direction == DIR_INCOMING)
        {
            /** IMPORTANT !!! Since we have hugely shift the nt value, we make sure to use a long enough integer. */
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ((graine >> 2 )  + ( nt << ((_span-1)*2)) ) & _mask; // previous kmer
                T reverse = revcomp (forward, _span);

                if (forward < reverse)
                {
                    if (this->contains (forward))
                    {
                        edges[idx++].set (forward, STRAND_FORWARD, source.kmer, source.strand, (Nucleotide)nt, DIR_INCOMING);
                    }
                }
                else
                {
                    if (this->contains (reverse))
                    {
                        edges[idx++].set (reverse, STRAND_REVCOMP, source.kmer, source.strand, (Nucleotide)nt, DIR_INCOMING);
                    }
                }
            }
        }

        else  {  throw system::Exception ("Unknown direction for getting edges");  }

        /** We update the size of the container according to the number of found items. */
        edges.setSize (idx);

        /** We return the number of found items. */
        return idx;
    }

    /**********************************************************************/
    size_t getNodes (const Node<T>& source, NodeSet<T>& nodes, Direction direction)
    {
        size_t idx = 0;

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        T graine = ((source.strand == STRAND_FORWARD) ?  source.kmer : revcomp (source.kmer, _span));

        if (direction == DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ( (graine << 2 )  + nt) & _mask;    // next kmer
                T reverse = revcomp (forward, _span);

                if (forward < reverse)  {  if (this->contains (forward))  {  nodes[idx++].set (forward, STRAND_FORWARD);  }  }
                else                    {  if (this->contains (reverse))  {  nodes[idx++].set (reverse, STRAND_REVCOMP);  }  }
            }
        }
        else if (direction == DIR_INCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ((graine >> 2 )  + ( nt <<  ((_span-1)*2)) ) & _mask; // previous kmer
                T reverse = revcomp (forward, _span);

                //std::cout << "span=" << _span << "  graine=" << graine << "  nt=" << (int)nt << "  forward=" << forward << "  reverse=" << reverse << std::endl;

                if (forward < reverse)  {  if (this->contains (forward))  {  nodes[idx++].set (forward, STRAND_FORWARD);  }  }
                else                    {  if (this->contains (reverse))  {  nodes[idx++].set (reverse, STRAND_REVCOMP);  }  }
            }
        }

        else  {  throw system::Exception ("Unknown direction for getting nodes");  }

        /** We update the size of the container according to the number of found items. */
        nodes.setSize (idx);

        /** We return the number of found items. */
        return idx;
    }

    /**********************************************************************/
    class NodeIteratorBasic : public INodeIterator<T>
    {
    public:

        NodeIteratorBasic (GraphBasic* graph, tools::dp::Iterator<T>* ref) : _ref(0), _current(graph)
        {
            setRef(ref);
        }

        ~NodeIteratorBasic ()  { setRef(0);   }

        /** \copydoc  Iterator::first */
        void first() { _ref->first(); }

        /** \copydoc  Iterator::next */
        void next()  { _ref->next(); }

        /** \copydoc  Iterator::isDone */
        bool isDone() { return _ref->isDone();  }

        /** \copydoc  Iterator::item */
        Node<T>& item ()  {  _current.kmer = _ref->item();  return _current; }

    private:
        tools::dp::Iterator<T>* _ref;
        void setRef (tools::dp::Iterator<T>* ref)  { SP_SETATTR(ref); }

        Node<T> _current;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
