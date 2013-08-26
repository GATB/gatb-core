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
    : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info(0)
{
    setProps (props);

    setInfo (new tools::misc::impl::Properties());

    setNodesIterable (new tools::collections::impl::IterableFile<T> (props->getStr (STR_KMER_SOLID)));

    setModel (new kmer::impl::Model<T> (props->getInt(STR_KMER_SIZE)));

    /** We compute the bloom size. */
    float NBITS_PER_KMER = log (16*props->getInt(STR_KMER_SIZE)*(lg2*lg2))/(lg2*lg2);
    size_t nbHash        = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t bloomSize =  (u_int64_t) (system::impl::System::file().getSize (props->getStr (STR_KMER_SOLID)) / sizeof (T) * NBITS_PER_KMER);

    setBloom (kmer::impl::BloomBuilder<T> (_nodesIterable->iterator(), bloomSize, nbHash).build ());

    if (props->get(STR_KMER_CFP) != 0)
    {
        setcFPset (new tools::collections::impl::ContainerSet<T> (
            new tools::collections::impl::IteratorFile<T> (props->getStr(STR_KMER_CFP))
        ));
    }

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
GraphBasic<T>::GraphBasic (
    tools::collections::Iterable<T>* solidKmers,
    tools::collections::Iterable<T>* cFPKmers,
    size_t kmerSize
)
    : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info(0)
{
    setInfo (new tools::misc::impl::Properties());

    setNodesIterable (solidKmers);

    setModel (new kmer::impl::Model<T> (kmerSize));

    /** We compute the bloom size. */
    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
    size_t nbHash        = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t bloomSize =  (u_int64_t) (solidKmers->getNbItems() * NBITS_PER_KMER);

    setBloom (kmer::impl::BloomBuilder<T> (_nodesIterable->iterator(), bloomSize, nbHash).build ());

    setcFPset (new tools::collections::impl::ContainerSet<T> (cFPKmers->iterator()));

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
GraphBasic<T>::GraphBasic (bank::IBank* bank, size_t kmerSize, size_t nks)
    : _nodesIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info(0)
{
    setInfo (new tools::misc::impl::Properties());
    _info->add (0, "graph");

    /** We set the model. */
    setModel (new kmer::impl::Model<T> (kmerSize));
    _mask = _model->getMask ();
    _span = _model->getSpan ();

    /** We create a DSK instance and execute it. */
    DSKAlgorithm<T> dsk (bank, kmerSize, nks);
    dsk.execute();
    getInfo()->add (1, dsk.getInfo());

    /** We create a debloom instance and execute it. */
    DebloomAlgorithm<T> debloom (dsk.getSolidKmers(), kmerSize);
    debloom.execute();
    getInfo()->add (1, debloom.getInfo());

    /** We set the solid and critical kmers iterables. */
    setNodesIterable (dsk.getSolidKmers());
    setcFPset        (new tools::collections::impl::ContainerSet<T> (debloom.getCriticalKmers()->iterator()));

    /** We compute the bloom size. */
    float NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
    size_t nbHash        = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t bloomSize =  (u_int64_t) (dsk.getSolidKmers()->getNbItems() * NBITS_PER_KMER);

    DEBUG (("nbSolid=%d  NBITS_PER_KMER=%f  bloomSize=%d \n", dsk.getSolidKmers()->getNbItems(), NBITS_PER_KMER, bloomSize));

    /** We build the bloom filter. */
    setBloom (kmer::impl::BloomBuilder<T> (dsk.getSolidKmers()->iterator(), bloomSize, nbHash).build ());
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
    setInfo          (0);
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
