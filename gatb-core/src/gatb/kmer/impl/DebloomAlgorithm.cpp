/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/math/NativeInt8.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Debloom: read solid kmers            ";
static const char* progressFormat2 = "Debloom: build extension             ";
static const char* progressFormat3 = "Debloom: finalization                ";

/********************************************************************************/

template <typename Item>
class BuildKmerExtension : public IteratorFunctor
{
public:
    void operator() (const Kmer<Item>& kmer)
    {
        /** We configure the neighbor kmers iterator for a given source kmer. */
        _itNeighbors.setSource (kmer.value);

        /** We iterate all 8 neighbors. */
        for (_itNeighbors.first(); !_itNeighbors.isDone(); _itNeighbors.next())
        {
            /** If the bloom contains the current neighbor, we add it to the debloom file. */
            if (_bloom->contains (*_itNeighbors))
            {
                _extendBag.insert (*_itNeighbors);
            }
        }
    }

    BuildKmerExtension (Model<Item>& model, Bloom<Item>* bloom, Bag<Item>* extendBag)
        : _bloom(bloom), _extendBag(extendBag, 50*1000, this->newSynchro()), _itNeighbors(model)  { }

private:
    Bloom<Item>*   _bloom;
    BagCache<Item> _extendBag;
    typename ModelAbstract<Item>::KmerNeighborIterator _itNeighbors;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
DebloomAlgorithm<ProductFactory,T>::DebloomAlgorithm (
    Product<ProductFactory>& product,
    Iterable<Kmer<T> >* solidIterable,
    size_t              kmerSize,
    BloomFactory::Kind  bloomKind,
    const std::string&  debloomUri,
    size_t              max_memory,
    IProperties*        options
)
    :  Algorithm("debloom", options), _product(product), _kmerSize(kmerSize), _bloomKind(bloomKind), _debloomUri("debloom"),
       _max_memory(max_memory),
       _solidIterable(0)
{
    setSolidIterable    (solidIterable);

    /** We get a collection for the cFP from the product. */
    setCriticalCollection (& product().template addCollection<T> ("debloom/cfp"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
DebloomAlgorithm<ProductFactory,T>::~DebloomAlgorithm ()
{
    setSolidIterable      (0);
    setCriticalCollection (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
void DebloomAlgorithm<ProductFactory,T>::execute ()
{
    Model<T> model (_kmerSize);

    /***************************************************/
    /** We create a bloom and insert solid kmers into. */
    /***************************************************/
    IProperties* bloomProps = new Properties();  LOCAL (bloomProps);
    Bloom<T>* bloom = createBloom (_solidIterable, bloomProps);
    bloom->use ();

    /*************************************************/
    /** We build the solid neighbors extension.      */
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "fill debloom file");

        /** We create an iterator over the solid kmers. */
        Iterator<Kmer<T> >* itKmers = createIterator<Kmer<T> > (
            _solidIterable->iterator(),
            _solidIterable->getNbItems(),
            progressFormat2
        );
        LOCAL (itKmers);

        /** We iterate the solid kmers and build the neighbors extension. */
        getDispatcher()->iterate (itKmers, BuildKmerExtension<T> (model, bloom, new BagFile<T>(_debloomUri)));
    }

    /** We save the bloom. */
    Collection<NativeInt8>* bloomCollection = & _product().template addCollection<NativeInt8> ("debloom/bloom");
    bloomCollection->insert ((NativeInt8*)bloom->getArray(), bloom->getSize());
    bloomCollection->addProperty ("properties", bloomProps->getXML());

    /** We get rid of the bloom. */
    bloom->forget ();

    /*************************************************************/
    /** We extract the solid kmers from the neighbors extension. */
    /*************************************************************/
    string inputUri  = _debloomUri;
    string outputUri = _debloomUri + "2";

    /** We need a hash that will hold solid kmers. */
    Hash16<T> partition (_max_memory);

    {
        TIME_INFO (getTimeInfo(), "finalize debloom file");

        Iterator<Kmer<T> >* itKmers = createIterator<Kmer<T> > (
            _solidIterable->iterator(),
            _solidIterable->getNbItems(),
            progressFormat3
        );
        LOCAL (itKmers);

        /** We iterate the solid kmers. */
        for (itKmers->first(); !itKmers->isDone(); itKmers->next())
        {
            /** We add the current kmer into the partition. */
            partition.insert (itKmers->item().value);

            /** We may have reach the maximum number of items in the partition. */
            if (partition.size() >= partition.getMaxNbItems())
            {
                /** We exclude the partition content from the critical false positive file. */
                end_debloom_partition (
                    partition,
                    new IteratorFile <T> (inputUri),
                    new BagFile      <T> (outputUri)
                );

                /** We swap the filenames. */
                std::swap (inputUri, outputUri);
            }
        }

        /** We finally write into the dedicated bag. */
        end_debloom_partition (
            partition,
            new IteratorFile <T> (inputUri),
            _criticalCollection->bag()
        );
    }

    /** We remove the two temporary files. */
    System::file().remove (inputUri);
    System::file().remove (outputUri);

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "critical kmers nb", "%ld", _criticalCollection->iterable()->getNbItems() );
    getInfo()->add (2, "critical kmers uri",       _criticalCollection->getFullId());
    getInfo()->add (1, getTimeInfo().getProperties("time"));
}

/*********************************************************************/

/** The following functor builds the critical false positive file by excluding
 *  kmers from a current cFP file if they are solid kmers.
 *  
 *  It is designed to get a kmer from the current cFP file and check
 *  whether it belongs to the solid kmers file (at least one part of it
 *  as a partition hash table).
 *
 *  In case the input cFP kmer doesn't belong to the solid kmers,
 *  we add it to the output cFP file.
 *
 *  This functor is supposed to be used concurrently by several thread. Therefore
 *  we have to protect the output cFP file against concurrent access, which is
 *  achieved by encapsulating the actual output file by a BagCache instance.
 */
template <typename Item>
class EndDebloomPartition : public IteratorFunctor
{
public:
    void operator() (const Item& extensionKmer)
    {
        /** If the extension kmer is not a solid kmer, then this is a critical false positive. */
        if (_solidKmers.contains (extensionKmer) == false)  {  _output.insert (extensionKmer);  }
    }

    EndDebloomPartition (Hash16<Item>& set, Bag<Item>* output)   : _solidKmers(set), _output (output, 10*1000, this->newSynchro())  { }

private:
    /** Reference of the solid kmers (as a hash table). */
    Hash16<Item>&   _solidKmers;

    /** Encapsulation of the output cFP file => ensure caching and thread-access protection. */
    BagCache<Item>  _output;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
void DebloomAlgorithm<ProductFactory,T>::end_debloom_partition (
    Hash16<T>&    partition,
    Iterator<T>*  debloomInput,
    Bag<T>*       debloomOutput
)
{
    LOCAL (debloomInput);
    LOCAL (debloomOutput);

    /** We dispatch the iteration of the current cFP bag. */
    getDispatcher()->iterate (debloomInput, EndDebloomPartition<T> (partition, debloomOutput));

    /** We make sure we output all the items. */
    debloomOutput->flush ();

    /** We clear the set. */
    partition.clear ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
Bloom<T>* DebloomAlgorithm<ProductFactory,T>::createBloom (
    tools::collections::Iterable<Kmer<T> >* solidIterable,
    tools::misc::IProperties* bloomProps
)
{
    TIME_INFO (getTimeInfo(), "create bloom from kmers");

    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

    u_int64_t solidKmersNb = solidIterable->getNbItems();

    u_int64_t estimatedBloomSize = (u_int64_t) (solidKmersNb * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator<Kmer<T> >* itKmers = createIterator<Kmer<T> > (
        solidIterable->iterator(),
        solidKmersNb,
        progressFormat1
    );
    LOCAL (itKmers);

    /** We use a bloom builder. */
    BloomBuilder<T> builder (itKmers, estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER), _bloomKind, getDispatcher()->getExecutionUnitsNumber());

    /** We instantiate the bloom object. */
    Bloom<T>* bloom = builder.build (bloomProps);
    bloomProps->add (1, "nbits per kmer", "%f", NBITS_PER_KMER);

    /** We return the created bloom filter. */
    return bloom;
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomAlgorithm <ProductFileFactory, NativeInt64>;
#ifdef INT128_FOUND
template class DebloomAlgorithm <ProductFileFactory, NativeInt128>;
#else
template class DebloomAlgorithm <ProductFileFactory, LargeInt<2> >;
#endif

template class DebloomAlgorithm <ProductFileFactory, LargeInt<3> >;
template class DebloomAlgorithm <ProductFileFactory, LargeInt<4> >;

/********************************************************************************/

template class DebloomAlgorithm <ProductHDF5Factory, NativeInt64>;
#ifdef INT128_FOUND
template class DebloomAlgorithm <ProductHDF5Factory, NativeInt128>;
#else
template class DebloomAlgorithm <ProductHDF5Factory, LargeInt<2> >;
#endif

template class DebloomAlgorithm <ProductHDF5Factory, LargeInt<3> >;
template class DebloomAlgorithm <ProductHDF5Factory, LargeInt<4> >;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
