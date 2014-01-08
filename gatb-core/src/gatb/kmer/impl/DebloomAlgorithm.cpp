/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

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

using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Debloom: read solid kmers              ";
static const char* progressFormat2 = "Debloom: build extension               ";
static const char* progressFormat3 = "Debloom: finalization                  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
DebloomAlgorithm<span>::DebloomAlgorithm (
    Product& product,
    Iterable<Count>*    solidIterable,
    size_t              kmerSize,
    size_t              max_memory,
    size_t              nb_cores,
    BloomFactory::Kind  bloomKind,
    const std::string&  debloomUri,
    IProperties*        options
)
    :  Algorithm("debloom", nb_cores, options), _product(product), _kmerSize(kmerSize), _bloomKind(bloomKind), _debloomUri("debloom"),
       _max_memory(max_memory),
       _solidIterable(0)
{
    /** We get a group for deblooming. */
    Group& group = _product().getGroup ("debloom");

    setSolidIterable    (solidIterable);

    /** We get a collection for the cFP from the product. */
    setCriticalCollection (& group.template getCollection<Type> ("cfp"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
DebloomAlgorithm<span>::~DebloomAlgorithm ()
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
template<size_t span>
void DebloomAlgorithm<span>::execute ()
{
    Model model (_kmerSize);

    /** We get a group for deblooming. */
    Group& group = _product().getGroup ("debloom");

    /***************************************************/
    /** We create a bloom and insert solid kmers into. */
    /***************************************************/
    IProperties* bloomProps = new Properties();  LOCAL (bloomProps);
    Bloom<Type>* bloom = createBloom (_solidIterable, bloomProps);
    bloom->use ();

    /*************************************************/
    /** We build the solid neighbors extension.      */
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "fill_debloom_file");

        /** We create an iterator over the solid kmers. */
        Iterator<Count>* itKmers = createIterator<Count> (
            _solidIterable->iterator(),
            _solidIterable->getNbItems(),
            progressFormat2
        );
        LOCAL (itKmers);

        /** We create a synchronized cache on the debloom output. This cache will be cloned by the dispatcher. */
        ThreadObject<BagCache<Type> > extendBag = BagCache<Type> (new BagFile<Type>(_debloomUri), 50*1000, System::thread().newSynchronizer());

        /** We iterate the solid kmers. */
        getDispatcher()->iterate (itKmers, [&] (const Count& kmer)
        {
            /** We iterate the neighbors of the current solid kmer. */
            model.iterateNeighbors (kmer.value, [&] (const Type& k)
            {
                if (bloom->contains (k))  {  extendBag().insert (k);  }
            });
        });

        /** We have to flush each bag cache used during iteration. */
        extendBag.foreach ([] (BagCache<Type>& bag)  { bag.flush(); });
    }

    /** We save the bloom. */
    Collection<NativeInt8>* bloomCollection = & group.template getCollection<NativeInt8> ("bloom");
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
    Hash16<Type> partition (_max_memory);

    {
        TIME_INFO (getTimeInfo(), "finalize_debloom_file");

        Iterator<Count>* itKmers = createIterator<Count> (
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
                    new IteratorFile <Type> (inputUri),
                    new BagFile      <Type> (outputUri)
                );

                /** We swap the filenames. */
                std::swap (inputUri, outputUri);
            }
        }

        /** We finally write into the dedicated bag. */
        end_debloom_partition (
            partition,
            new IteratorFile <Type> (inputUri),
            _criticalCollection->bag()
        );
    }

    /** We remove the two temporary files. */
    System::file().remove (inputUri);
    System::file().remove (outputUri);

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, bloomProps);
    getInfo()->add (2, "critical_kmers_nb", "%ld", _criticalCollection->iterable()->getNbItems() );

    getInfo()->add (1, getTimeInfo().getProperties("time"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void DebloomAlgorithm<span>::end_debloom_partition (
    Hash16<Type>&    partition,
    Iterator<Type>*  debloomInput,
    Bag<Type>*       debloomOutput
)
{
    LOCAL (debloomInput);
    LOCAL (debloomOutput);

    /** We create a synchronized cache on the debloom output. This cache will be cloned by the dispatcher. */
    ThreadObject<BagCache<Type> > finalizeBag = BagCache<Type> (debloomOutput, 8*1024, System::thread().newSynchronizer());

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
    getDispatcher()->iterate (debloomInput, [&] (const Type& extensionKmer)
    {
        if (partition.contains (extensionKmer) == false)  {  finalizeBag().insert (extensionKmer);  }
    });

    /** We have to flush each bag cache used during iteration. */
    finalizeBag.foreach ([] (BagCache<Type>& bag)  { bag.flush(); });

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
template<size_t span>
Bloom <typename Kmer<span>::Type>* DebloomAlgorithm<span>::createBloom (
    tools::collections::Iterable<Count>* solidIterable,
    tools::misc::IProperties* bloomProps
)
{
    TIME_INFO (getTimeInfo(), "create_bloom_from_kmers");

    /** We get the number of solid kmers. */
    u_int64_t solidKmersNb = solidIterable->getNbItems();

    double lg2 = log(2);
    float     NBITS_PER_KMER     = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    size_t    nbHash             = (int)floorf (0.7*NBITS_PER_KMER);
    u_int64_t estimatedBloomSize = (u_int64_t) (solidIterable->getNbItems() * NBITS_PER_KMER);

    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator <Count>* itKmers = createIterator<Count> (
        solidIterable->iterator(),
        solidKmersNb,
        progressFormat1
    );
    LOCAL (itKmers);

    /** We use a bloom builder. */
    BloomBuilder<span> builder (estimatedBloomSize, nbHash, _bloomKind, getDispatcher()->getExecutionUnitsNumber());

    /** We instantiate the bloom object. */
    Bloom<Type>* bloom = builder.build (itKmers, bloomProps);
    bloomProps->add (0, "nbits_per_kmer", "%f", NBITS_PER_KMER);

    /** We return the created bloom filter. */
    return bloom;
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomAlgorithm <32>;
template class DebloomAlgorithm <64>;
template class DebloomAlgorithm <96>;
template class DebloomAlgorithm <128>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
