/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/debruijn/impl/GraphBasic.hpp>

#include <gatb/kmer/impl/DSKAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <gatb/tools/math/NativeInt8.hpp>

#include <gatb/system/impl/System.hpp>

using namespace std;
using namespace gatb::core::tools::math;

using namespace gatb::core::kmer;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::kmer::impl;

using namespace gatb::core::bank::impl;

using namespace gatb::core::system::impl;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename T>
GraphBasic<T>::GraphBasic (tools::misc::IProperties* options)
    : _state(INIT), _solidKmersIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info("graph"), _product(0)
{
    /** We keep a reference on some provided parameters. */
    setProps (options);
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
    tools::collections::Iterable<Kmer<T> >* solidKmers,
    tools::collections::Iterable<T>*        cFPKmers,
    size_t kmerSize
)
    : _state(INIT), _solidKmersIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info("graph"), _product(0)
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
GraphBasic<T>::GraphBasic (const std::string& uri)
    : _state(INIT), _solidKmersIterable(0), _bloom(0), _cFPset(0), _model(0), _props(0), _info("graph"), _product(0)
{
    /** We create a product instance. */
    setProduct (ProductFactoryLocal::createProduct (uri, false, false));

    /** We set the iterable for the solid kmers. */
    tools::collections::Collection<Kmer<T> >* solidCollection = & getProduct("dsk").template getCollection<Kmer<T> > ("solid");
    setSolidKmersIterable (solidCollection);

    /** We retrieve the properties of the solid collection. */
    string solidPropsString = solidCollection->getProperty("properties");
    stringstream ss (solidPropsString);

    /** We build a Properties instance from the retrieved XML string. */
    Properties* solidProps = new Properties ();  LOCAL (solidProps);
    solidProps->readXML (ss);

    /** We get the kmer size from the properties. */
    IProperty* prop = solidProps->get ("kmer_size");

    if (prop == 0)  { throw "UNABLE TO GET PROPERTIES OF THE SOLID KMER DATASET"; }
    size_t kmerSize = solidProps->getInt ("kmer_size");

    /** We create the kmer model. */
    setModel (new kmer::impl::Model<T> (kmerSize));

    /** We compute parameters for the Bloom filter. */
    double lg2 = log(2);
    float     NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
    size_t    nbHash         = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t bloomSize      = (u_int64_t) (solidCollection->getNbItems() * NBITS_PER_KMER);

    tools::collections::Iterable<NativeInt8>* bloomIterable = & getProduct("debloom").template getCollection<NativeInt8> ("bloom");

    tools::collections::impl::Bloom<T>* bloom =  tools::collections::impl::BloomFactory::singleton().createBloom<T> (
        tools::collections::impl::BloomFactory::CacheCoherent,
        bloomSize,
        nbHash
    );

    setBloom (bloom);
    bloomIterable->getItems ((NativeInt8*&)_bloom->getArray());

    tools::collections::Iterable<T>* cFPKmers = & getProduct("debloom").template getCollection<T > ("cfp");
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
GraphBasic<T>::~GraphBasic ()
{
    setProps              (0);
    setSolidKmersIterable (0);
    setcFPset             (0);
    setModel              (0);
    setBloom              (0);
    setProduct            (0);

    _state = INIT;
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
void GraphBasic<T>::build (bank::IBank* bank)
{
    LOCAL (bank);

    if (getState() == INIT)
    {
        _kmerSize = _props->get(STR_KMER_SIZE)  ? _props->getInt(STR_KMER_SIZE) : 27;
        _nks      = _props->get(STR_NKS)        ? _props->getInt(STR_NKS)       : 3;

        string output   = _props->get(STR_URI_OUTPUT) ?
                _props->getStr(STR_URI_OUTPUT)   :
                system::impl::System::file().getBaseName (bank->getId());

        string binaryBankUri = System::file().getCurrentDirectory() + "/bank.bin";

        /************************************************************/
        /*                       Product creation                   */
        /************************************************************/
        setProduct (ProductFactoryLocal::createProduct (output, true, false));

        /************************************************************/
        /*                         Bank conversion                  */
        /************************************************************/
        /** We create the binary bank. */
        BankConverterAlgorithm converter (bank, _kmerSize, binaryBankUri);
        executeAlgorithm (converter);

        /************************************************************/
        /*                         DSK                              */
        /************************************************************/
        /** We create a DSK instance and execute it. */
        DSKAlgorithm<ProductFactoryLocal, T> dsk (
            *_product,
            converter.getResult(),
            _kmerSize,
            _nks,
            _props->get(STR_MAX_MEMORY) ? _props->getInt(STR_MAX_MEMORY) : 1000,
            _props->get(STR_MAX_DISK)   ? _props->getInt(STR_MAX_DISK)   : 0,
            _props->get(STR_NB_CORES)   ? _props->getInt(STR_NB_CORES)   : 0
        );
        executeAlgorithm (dsk);

        /************************************************************/
        /*                         Debloom                          */
        /************************************************************/
        /** We create a debloom instance and execute it. */
        DebloomAlgorithm<ProductFactoryLocal, T> debloom (
            *_product,
            dsk.getSolidKmers(),
            _kmerSize,
            _props->get(STR_MAX_MEMORY) ? _props->getInt(STR_MAX_MEMORY) : 1000,
            _props->get(STR_NB_CORES)   ? _props->getInt(STR_NB_CORES)   : 0
        );
        executeAlgorithm (debloom);

        /** We retrieve the bloom raw data. */
        Iterable<NativeInt8>* bloomArray = & getProduct("debloom").template getCollection<NativeInt8> ("bloom");

        /** We configure the graph from the result of DSK and debloom. */
        configure (_kmerSize, dsk.getSolidKmers(), debloom.getCriticalKmers(), bloomArray);

        /** We can get rid of the binary bank. */
        System::file().remove (binaryBankUri);

        /** We add metadata to some collections. */
        dsk.getSolidKmers()->addProperty ("properties", dsk.getInfo()->getXML());
        debloom.getCriticalKmers()->addProperty ("properties", debloom.getInfo()->getXML());

        /** We update the state of the graph. */
        setState (READY);
    }
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
    tools::collections::Iterable<Kmer<T> >*   solidIterable,
    tools::collections::Iterable<T>*          cFPKmers,
    tools::collections::Iterable<tools::math::NativeInt8>* bloomArray
)
{
    /** We create the kmer model. */
    setModel (new kmer::impl::Model<T> (kmerSize));

    /** We set the iterable for the solid kmers. */
    setSolidKmersIterable (solidIterable);

    /** We compute parameters for the Bloom filter. */
    double lg2 = log(2);
    float     NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
    size_t    nbHash         = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t bloomSize      = (u_int64_t) (solidIterable->getNbItems() * NBITS_PER_KMER);

    /** We need a bloom builder for creating the Bloom filter. */
    kmer::impl::BloomBuilder<T> builder (bloomSize, nbHash);

    /** We check whether we have a stored Bloom filter. */
    if (bloomArray == 0)
    {
        /** We create Bloom filter and fill it with the solid kmers. */
        setBloom (builder.build (solidIterable->iterator()));
    }
    else
    {
        /** We load the Bloom filter from the specific dataset. */
        setBloom (builder.load (bloomArray));
    }

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
void GraphBasic<T>::executeAlgorithm (gatb::core::tools::misc::impl::Algorithm& algorithm)
{
    string bargraph = _props->get(STR_VERBOSE) ?  "2" : "0";

    algorithm.getInput()->add (0, STR_PROGRESS_BAR, bargraph);

    algorithm.execute();

    getInfo()->add (1, algorithm.getInfo());
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
    T graine = ((source.strand == STRAND_FORWARD) ?  source.kmer.value : revcomp (source.kmer.value, _span));

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
                    edges[idx++].set (source.getGraph(), source.kmer.value, source.strand, forward, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                }
            }
            else
            {
                if (this->contains (reverse))
                {
                    edges[idx++].set (source.getGraph(), source.kmer.value, source.strand, reverse, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
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
                    edges[idx++].set (source.getGraph(), forward, STRAND_FORWARD, source.kmer.value, source.strand, (Nucleotide)nt, DIR_INCOMING);
                }
            }
            else
            {
                if (this->contains (reverse))
                {
                    edges[idx++].set (source.getGraph(), reverse, STRAND_REVCOMP, source.kmer.value, source.strand, (Nucleotide)nt, DIR_INCOMING);
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
    T graine = ((source.strand == STRAND_FORWARD) ?  source.kmer.value : revcomp (source.kmer.value, _span));

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

template class GraphBasic <NativeInt64>;
#ifdef INT128_FOUND
template class GraphBasic <NativeInt128>;
#else
template class GraphBasic <LargeInt<2> >;
#endif

template class GraphBasic <LargeInt<3> >;
template class GraphBasic <LargeInt<4> >;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
