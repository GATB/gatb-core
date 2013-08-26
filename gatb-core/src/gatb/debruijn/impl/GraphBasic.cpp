/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/debruijn/impl/GraphBasic.hpp>

#include <gatb/kmer/impl/DSKAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/tools/misc/impl/Property.hpp>

using namespace std;
using namespace gatb::core::tools::math;

using namespace gatb::core::kmer::impl;

#define DEBUG(a)  //printf a

bool dbg = false;
/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

static const double lg2 = 0.69314718055994530941723212145817656808;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
GraphBasic<T>::GraphBasic (tools::misc::IProperties* props)
    : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info("graph")
{
    setProps (props);

    assert (props->get(STR_KMER_SIZE)  != 0);
    assert (props->get(STR_KMER_SOLID) != 0);
    assert (props->get(STR_KMER_CFP)   != 0);

    configure (
        props->getInt(STR_KMER_SIZE),
        new tools::collections::impl::IterableFile<T> (props->getStr (STR_KMER_SOLID)),
        new tools::collections::impl::IterableFile<T> (props->getStr(STR_KMER_CFP))
    );
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
GraphBasic<T>::GraphBasic (
    tools::collections::Iterable<T>* solidKmers,
    tools::collections::Iterable<T>* cFPKmers,
    size_t kmerSize
)
    : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info("graph")
{

    configure (kmerSize, solidKmers, cFPKmers);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
GraphBasic<T>::GraphBasic (bank::IBank* bank, size_t kmerSize, size_t nks)
    : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info("graph")
{
    /** We create a DSK instance and execute it. */
    DSKAlgorithm<T> dsk (bank, kmerSize, nks);
    dsk.getInput()->add (0, STR_PROGRESS_BAR, "2");
    dsk.execute();
    getInfo()->add (1, dsk.getInfo());

    /** We create a debloom instance and execute it. */
    DebloomAlgorithm<T> debloom (dsk.getSolidKmers(), kmerSize);
    debloom.getInput()->add (0, STR_PROGRESS_BAR, "2");
    debloom.execute();
    getInfo()->add (1, debloom.getInfo());

    /** We configure the graph from the result of DSK and debloom. */
    configure (kmerSize, dsk.getSolidKmers(), debloom.getCriticalKmers());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
GraphBasic<T>::~GraphBasic ()
{
    setProps         (0);
    setNodesIterable (0);
    setcFPset        (0);
    setModel         (0);
    setBloom         (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
void GraphBasic<T>::configure (
    size_t kmerSize,
    tools::collections::Iterable<T>* solidIterable,
    tools::collections::Iterable<T>* cFPKmers
)
{
    /** We create the kmer model. */
    setModel (new kmer::impl::Model<T> (kmerSize));

    /** We set the iterable for the solid kmers. */
    setNodesIterable (solidIterable);

    /** We compute parameters for the Bloom filter. */
    double lg2 = log(2);
    float     NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
    size_t    nbHash         = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t bloomSize      = (u_int64_t) (solidIterable->getNbItems() * NBITS_PER_KMER);

    /** We create Bloom filter and fill it with the solid kmers. */
    setBloom (kmer::impl::BloomBuilder<T> (solidIterable->iterator(), bloomSize, nbHash).build ());

    /** We build the set of critical false positive kmers. */
    setcFPset (new tools::collections::impl::ContainerSet<T> (cFPKmers->iterator()));

    /** Some shortcuts attributes. */
    _mask = _model->getMask ();
    _span = _model->getSpan ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
size_t GraphBasic<T>::getEdges (const Node<T>& source, EdgeSet<T>& edges, Direction direction)
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
                    edges[idx++].set (source.getGraph(), source.kmer, source.strand, forward, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                }
            }
            else
            {
                if (this->contains (reverse))
                {
                    edges[idx++].set (source.getGraph(), source.kmer, source.strand, reverse, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
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
                    edges[idx++].set (source.getGraph(), forward, STRAND_FORWARD, source.kmer, source.strand, (Nucleotide)nt, DIR_INCOMING);
                }
            }
            else
            {
                if (this->contains (reverse))
                {
                    edges[idx++].set (source.getGraph(), reverse, STRAND_REVCOMP, source.kmer, source.strand, (Nucleotide)nt, DIR_INCOMING);
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

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
size_t GraphBasic<T>::getNodes (const Node<T>& source, NodeSet<T>& nodes, Direction direction)
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

            if (forward < reverse)  {  if (this->contains (forward))  {  nodes[idx++].set (source.getGraph(), forward, STRAND_FORWARD);  }  }
            else                    {  if (this->contains (reverse))  {  nodes[idx++].set (source.getGraph(), reverse, STRAND_REVCOMP);  }  }
        }
    }
    else if (direction == DIR_INCOMING)
    {
        for (u_int64_t nt=0; nt<4; nt++)
        {
            T forward = ((graine >> 2 )  + ( nt <<  ((_span-1)*2)) ) & _mask; // previous kmer
            T reverse = revcomp (forward, _span);

            if (forward < reverse)  {  if (this->contains (forward))  {  nodes[idx++].set (source.getGraph(), forward, STRAND_FORWARD);  }  }
            else                    {  if (this->contains (reverse))  {  nodes[idx++].set (source.getGraph(), reverse, STRAND_REVCOMP);  }  }
        }
    }

    else  {  throw system::Exception ("Unknown direction for getting nodes");  }

    /** We update the size of the container according to the number of found items. */
    nodes.setSize (idx);

    /** We return the number of found items. */
    return idx;
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class GraphBasic <gatb::core::tools::math::NativeInt64>;
#ifdef INT128_FOUND
template class GraphBasic <gatb::core::tools::math::NativeInt128>;
#else
template class GraphBasic <gatb::core::tools::math::LargeInt<2> >;
#endif

template class GraphBasic <gatb::core::tools::math::LargeInt<3> >;
template class GraphBasic <gatb::core::tools::math::LargeInt<4> >;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
