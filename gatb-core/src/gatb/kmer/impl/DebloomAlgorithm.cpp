/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <gatb/kmer/impl/DebloomAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.pri>

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

#include <gatb/debruijn/impl/ContainerNode.hpp>

#include <gatb/tools/storage/impl/StorageTools.hpp>

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
static const char* progressFormat4 = "Debloom: cascading                     ";

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
    Storage& storage,
    Iterable<Count>*    solidIterable,
    size_t              kmerSize,
    size_t              max_memory,
    size_t              nb_cores,
    BloomFactory::Kind  bloomKind,
    DebloomKind         cascadingKind,
    const std::string&  debloomUri,
    IProperties*        options
)
    :  Algorithm("debloom", nb_cores, options), _storage(storage), _group(storage().getGroup ("debloom")),
       _kmerSize(kmerSize),
       _bloomKind(bloomKind), _cascadingKind(cascadingKind),
       _debloomUri("debloom"),
       _max_memory(max_memory),
       _solidIterable(0), _criticalCollection(0), _container(0)
{
    /** We get a group for deblooming. */
    Group& group = _storage().getGroup ("debloom");

    setSolidIterable    (solidIterable);

    /** We get a collection for the cFP from the storage. */
    setCriticalCollection (& group.getCollection<Type> ("cfp"));

    /** We set the max memory to a default project value if not set. */
    if (_max_memory == 0)  {  _max_memory = System::info().getMemoryProject(); }
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
DebloomAlgorithm<span>::DebloomAlgorithm (tools::storage::impl::Storage& storage)
:  Algorithm("debloom", 0, 0), _storage(storage), _group(storage().getGroup ("debloom")),
   _kmerSize(0),
   _debloomUri("debloom"),
   _max_memory(0),
   _solidIterable(0), _criticalCollection(0), _container(0)
{
    /** We get a collection for the cFP from the storage. */
    setCriticalCollection (& _group.getCollection<Type> ("cfp"));

    /** We retrieve the cascading kind from the storage. */
    parse (_criticalCollection->getProperty("kind"), _cascadingKind);

    loadDebloomStructures (_storage);
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
    setDebloomStructures  (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
#ifndef WITH_LAMBDA_EXPRESSION
template<typename Model, typename Count, typename Type>
struct FunctorKmersExtension
{
    struct FunctorNeighbors
    {
        Bloom<Type>* bloom; ThreadObject<BagCache<Type> >& extendBag;
        FunctorNeighbors (Bloom<Type>* bloom, ThreadObject<BagCache<Type> >& extendBag)
            : bloom(bloom), extendBag(extendBag) {}
        void operator() (const Type& k)  const {  if (bloom->contains (k))  {  extendBag().insert (k);  }  }
    } functorNeighbors;

    Model& model;
    FunctorKmersExtension (Model& model, Bloom<Type>* bloom, ThreadObject<BagCache<Type> >& extendBag)
        : model(model), functorNeighbors(bloom,extendBag) {}
    void operator() (const Count& kmer) const
    {
        /** We iterate the neighbors of the current solid kmer. */
        model.iterateNeighbors (kmer.value, functorNeighbors);
    }
};
#endif
/*********************************************************************/

template<size_t span>
void DebloomAlgorithm<span>::execute ()
{
    Model model (_kmerSize);

    Collection<Type>* criticalCollection = new CollectionFile<Type> ("cfp");
    LOCAL (criticalCollection);

    /***************************************************/
    /** We create a bloom and insert solid kmers into. */
    /***************************************************/
    IProperties* bloomProps = new Properties();  LOCAL (bloomProps);
    u_int64_t totalSizeBloom = 0;

    Bloom<Type>* bloom = createBloom (_solidIterable, bloomProps, totalSizeBloom);
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

#ifdef WITH_LAMBDA_EXPRESSION
        auto functorKmers = [&] (const Count& kmer)
        {
            /** We iterate the neighbors of the current solid kmer. */
            model.iterateNeighbors (kmer.value, [&] (const Type& k)
            {
                if (bloom->contains (k))  {  extendBag().insert (k);  }
            });
        };
#else
        FunctorKmersExtension<Model,Count,Type> functorKmers (model, bloom, extendBag);
#endif
        /** We iterate the solid kmers. */
        getDispatcher()->iterate (itKmers, functorKmers);

        /** We have to flush each bag cache used during iteration. */
        for (size_t i=0; i<extendBag.size(); i++)  {  extendBag[i].flush();  }
    }

    /** We save the bloom. */
    StorageTools::singleton().saveBloom<Type> (_group, "bloom", bloom);

    /** We get rid of the bloom. */
    bloom->forget ();

    /*************************************************************/
    /** We extract the solid kmers from the neighbors extension. */
    /*************************************************************/
    string inputUri  = _debloomUri;
    string outputUri = _debloomUri + "2";

    /** We need a hash table that will hold solid kmers. */
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
            criticalCollection->bag()
        );
    }

    /*************************************************************/
    /** We build the final cFP container.                        */
    /*************************************************************/
    IProperties* cfpProps = new Properties();  LOCAL (cfpProps);
    u_int64_t totalSizeCFP = 0;

    {
        TIME_INFO (getTimeInfo(), "cascading");

        /** We build the container node. */
        createCFP (criticalCollection, cfpProps, totalSizeCFP);
    }

    /** We remove the two temporary files. */
    System::file().remove (inputUri);
    System::file().remove (outputUri);
    criticalCollection->remove ();

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "bloom_kind",    "%s",  BloomFactory::toString(_bloomKind));
    getInfo()->add (2, "debloom_kind",  "%s", toString(_cascadingKind).c_str());
    getInfo()->add (2, "size",           "%ld", totalSizeBloom + totalSizeCFP);
    getInfo()->add (2, "nbits_per_kmer", "%f",  (float)(totalSizeBloom + totalSizeCFP) / (float)_solidIterable->getNbItems());
    getInfo()->add (2, "bloom");
    getInfo()->add (3, bloomProps);
    getInfo()->add (2, cfpProps);

    getInfo()->add (1, getTimeInfo().getProperties("time"));

    /** We save (as metadata) information into the critical collection. */
    _criticalCollection->addProperty ("properties", getInfo()->getXML());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
#ifndef WITH_LAMBDA_EXPRESSION
template<typename Type>
struct FunctorKmersFinalize
{
    Hash16<Type>& partition; ThreadObject<BagCache<Type> >& finalizeBag;
    FunctorKmersFinalize (Hash16<Type>& partition, ThreadObject<BagCache<Type> >& finalizeBag)
        : partition(partition), finalizeBag(finalizeBag) {}
    void operator() (const Type& extensionKmer)
    {
        if (partition.contains (extensionKmer) == false)  {  finalizeBag().insert (extensionKmer);  }
    }
};
#endif
/*********************************************************************/

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
#ifdef WITH_LAMBDA_EXPRESSION
    auto functorKmers = [&] (const Type& extensionKmer)
    {
        if (partition.contains (extensionKmer) == false)  {  finalizeBag().insert (extensionKmer);  }
    };
#else
    FunctorKmersFinalize<Type> functorKmers (partition, finalizeBag);
#endif

    getDispatcher()->iterate (debloomInput, functorKmers);

    /** We have to flush each bag cache used during iteration. */
    for (size_t i=0; i<finalizeBag.size(); i++) { finalizeBag[i].flush(); }

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
    tools::misc::IProperties* bloomProps,
    u_int64_t& totalSizeBloom
)
{
    TIME_INFO (getTimeInfo(), "create_bloom_from_kmers");

    /** We get the number of solid kmers. */
    u_int64_t solidKmersNb = solidIterable->getNbItems();

    float     NBITS_PER_KMER     = getNbBitsPerKmer();
    u_int64_t estimatedBloomSize = (u_int64_t) (solidKmersNb * NBITS_PER_KMER);
    size_t    nbHash             = (int)floorf (0.7*NBITS_PER_KMER);

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

    totalSizeBloom = bloom->getBitSize();

    /** We return the created bloom filter. */
    return bloom;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/

#ifndef WITH_LAMBDA_EXPRESSION

/*********************************************************************/
template<typename Type>
struct Functor0
{
    Bloom<Type>* bloom2;
    Functor0 (Bloom<Type>* bloom2) : bloom2(bloom2) {}
    void operator() (const Type& t)  {  bloom2->insert (t); };
};
/*********************************************************************/
template<typename Count, typename Type>
struct Functor1
{
    Iterable<Count>* solidIterable;  Bloom<Type>* bloom2;  Bloom<Type>* bloom3;  ThreadObject<BagCache<Type> >& T2Cache;
    Functor1 (Iterable<Count>* solidIterable,  Bloom<Type>* bloom2, Bloom<Type>* bloom3, ThreadObject<BagCache<Type> >& T2Cache)
        : solidIterable(solidIterable), bloom2(bloom2), bloom3(bloom3), T2Cache(T2Cache) {}
    void operator() (const Count& t)
    {
        if (bloom2->contains(t.value))
        {
            T2Cache().insert (t.value);
            bloom3->insert (t.value);
        }
    }
};
/*********************************************************************/
template<typename Type>
struct Functor2
{
    Bloom<Type>* bloom3;  Bloom<Type>* bloom4;
    Functor2 (Bloom<Type>* bloom3, Bloom<Type>* bloom4) : bloom3(bloom3), bloom4(bloom4) {}
    void operator() (const Type& t) { if (bloom3->contains(t))  {  bloom4->insert (t); }  }
};
/*********************************************************************/

#endif  // ! WITH_LAMBDA_EXPRESSION

/*********************************************************************/

template<size_t span>
void DebloomAlgorithm<span>::createCFP (
    Collection<Type>*  criticalCollection,
    IProperties* props,
    u_int64_t& totalSize_bits
)
{
    props->add (0, "cfp", "");

    /** We define an iterator for 4 tasks. */
    size_t nbTasks = 4;
    Iterator<int>* itTask = createIterator<int> (new Range<int>::Iterator (1,nbTasks), nbTasks, progressFormat4);
    LOCAL (itTask);
    itTask->first ();

    /** We may have to change the cascading kind if we have no false positives. */
    int64_t nbFP = criticalCollection->getNbItems();
    if (nbFP == 0)  {  _cascadingKind = DEBLOOM_ORIGINAL; }

    /*************************************************************/
    /** We code the critical FP according to the cascading type. */
    /*************************************************************/
    switch (_cascadingKind)
    {
        case DEBLOOM_CASCADING:
        {
            float     nbBitsPerKmer = getNbBitsPerKmer();
            u_int64_t nbKmers       = _solidIterable->getNbItems();

            int64_t estimated_T2_size = std::max ((int)ceilf(nbKmers * (double)powf((double)0.62, (double)nbBitsPerKmer)), 1);
            int64_t estimated_T3_size = std::max ((int)ceilf(nbFP    * (double)powf((double)0.62, (double)nbBitsPerKmer)), 1);

            tools::collections::impl::Bloom<Type>* bloom2 = tools::collections::impl::BloomFactory::singleton().createBloom<Type> (
                _bloomKind, (u_int64_t)(nbFP * nbBitsPerKmer), (int)floorf(0.7*nbBitsPerKmer)
            );
            LOCAL (bloom2);

            tools::collections::impl::Bloom<Type>* bloom3 = tools::collections::impl::BloomFactory::singleton().createBloom<Type> (
                _bloomKind, (u_int64_t)(estimated_T2_size * nbBitsPerKmer), (int)floorf(0.7*nbBitsPerKmer)
            );
            LOCAL (bloom3);

            tools::collections::impl::Bloom<Type>* bloom4 = tools::collections::impl::BloomFactory::singleton().createBloom<Type> (
                _bloomKind, (u_int64_t)(estimated_T3_size * nbBitsPerKmer), (int)floorf(0.7*nbBitsPerKmer)
            );
            LOCAL (bloom4);

            // **** Insert the false positives in B2 ****
#ifdef WITH_LAMBDA_EXPRESSION
            auto functor0 = [&] (const Type& t)  {  bloom2->insert (t); };
#else
            Functor0<Type> functor0 (bloom2);
#endif
            getDispatcher()->iterate (criticalCollection->iterator(), functor0);
            itTask->next();

            //  **** Insert false positives in B3 and write T2
            const char* T2name = "t2_kmers";
            Collection<Type>* T2File = new CollectionFile<Type> (T2name);  LOCAL (T2File);

            /** We need to protect the T2File against concurrent accesses, we use a ThreadObject for this. */
            ThreadObject<BagCache<Type> > T2Cache = BagCache<Type> (T2File, 8*1024, System::thread().newSynchronizer());

#ifdef WITH_LAMBDA_EXPRESSION
            auto functor1 = [&] (const Count& t)
            {
                if (bloom2->contains(t.value))
                {
                    T2Cache().insert (t.value);
                    bloom3->insert (t.value);
                }
            };
#else
            Functor1<Count,Type> functor1 (_solidIterable, bloom2, bloom3, T2Cache);
#endif
            getDispatcher()->iterate (_solidIterable->iterator(), functor1);

            /** We have to flush each bag cache used during iteration. */
            for (size_t i=0; i<T2Cache.size(); i++)  { T2Cache[i].flush(); }
            itTask->next();

            // **** Insert false positives in B4 (we could write T3, but it's not necessary)
#ifdef WITH_LAMBDA_EXPRESSION
            auto functor2 = [&] (const Type& t)  {  if (bloom3->contains(t))  {  bloom4->insert (t); }  };
#else
            Functor2<Type> functor2 (bloom3,bloom4);
#endif
            getDispatcher()->iterate (criticalCollection->iterator(), functor2);
            itTask->next();

            /** We build the final cfp set. */
            vector<Type> cfpItems;
            Iterator<Type>* T2It = T2File->iterator(); LOCAL(T2It);
            for (T2It->first(); !T2It->isDone(); T2It->next())
            {
                if (bloom4->contains(T2It->item()))  {  cfpItems.push_back (T2It->item());  }
            }
            std::sort (cfpItems.begin(), cfpItems.end());
            _criticalCollection->insert (cfpItems.data(), cfpItems.size());
            itTask->next();
            itTask->isDone(); // force to finish progress dump

            /** We save the final cFP container into the storage. */
            StorageTools::singleton().saveBloom<Type>      (_storage().getGroup ("debloom"), "bloom2", bloom2);
            StorageTools::singleton().saveBloom<Type>      (_storage().getGroup ("debloom"), "bloom3", bloom3);
            StorageTools::singleton().saveBloom<Type>      (_storage().getGroup ("debloom"), "bloom4", bloom4);

            totalSize_bits = bloom2->getBitSize() + bloom3->getBitSize() + bloom4->getBitSize() + 8*cfpItems.size()*sizeof(Type);

            /** Some statistics. */
            props->add (1, "size",   "%ld", totalSize_bits);
            props->add (2, "bloom2", "%ld", bloom2->getBitSize());
            props->add (2, "bloom3", "%ld", bloom3->getBitSize());
            props->add (2, "bloom4", "%ld", bloom4->getBitSize());
            props->add (2, "cfp",    "%ld", 8*cfpItems.size()*sizeof(Type));

            /** We remove the two temporary files. */
            System::file().remove (T2name);

            break;
        }

        case DEBLOOM_ORIGINAL:
        case DEBLOOM_DEFAULT:
        default:
        {
            /** We save the final cFP container into the storage. */
            tools::dp::Iterator<Type>* it = criticalCollection->iterator();  LOCAL (it);
            for (it->first(); !it->isDone(); it->next())  {  this->_criticalCollection->insert (it->item());  }
            _criticalCollection->flush ();

            totalSize_bits = 8*criticalCollection->getNbItems()*sizeof(Type);

            props->add (1, "size", "%ld", totalSize_bits);

            break;
        }
    }

    /** We save the cfp kind in the "cfp" storage collection. */
    _criticalCollection->addProperty ("kind", toString(_cascadingKind));
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
float DebloomAlgorithm<span>::getNbBitsPerKmer () const
{
    static double lg2 = log(2);
    float nbitsPerKmer = 0;

    if (_kmerSize > 128 && _cascadingKind==DEBLOOM_CASCADING)  {  throw Exception ("kmer size %d too big for cascading bloom filters", _kmerSize); }

    switch (_cascadingKind)
    {
    case DEBLOOM_CASCADING:
        nbitsPerKmer = rvalues[_kmerSize][1];
        break;

    case DEBLOOM_ORIGINAL:
    case DEBLOOM_DEFAULT:
    default:
        nbitsPerKmer = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
        break;
    }
    return nbitsPerKmer;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : used to be named loadContainer but I renamed it for clarity, also conficting name with actual loadContainer for a container
*********************************************************************/
template<size_t span>
void DebloomAlgorithm<span>::loadDebloomStructures (tools::storage::impl::Storage& storage)
{
    switch (_cascadingKind)
    {
        case DEBLOOM_ORIGINAL:
        case DEBLOOM_DEFAULT:
        default:
        {
            Bloom<Type>*       bloom    = StorageTools::singleton().loadBloom<Type>     (_group, "bloom");
            Container<Type>*   cFP      = StorageTools::singleton().loadContainer<Type> (_group, "cfp");

            /** We build the set of critical false positive kmers. */
            setDebloomStructures (new debruijn::impl::ContainerNode<Type> (bloom, cFP));

            break;
        }

        case DEBLOOM_CASCADING:
        {
            Bloom<Type>*     bloom   = StorageTools::singleton().loadBloom<Type>     (_group, "bloom");
            Bloom<Type>*     bloom2  = StorageTools::singleton().loadBloom<Type>     (_group, "bloom2");
            Bloom<Type>*     bloom3  = StorageTools::singleton().loadBloom<Type>     (_group, "bloom3");
            Bloom<Type>*     bloom4  = StorageTools::singleton().loadBloom<Type>     (_group, "bloom4");
            Container<Type>* cFP     = StorageTools::singleton().loadContainer<Type> (_group, "cfp");

            /** We build the set of critical false positive kmers. */
            setDebloomStructures (new debruijn::impl::ContainerNodeCascading<Type> (bloom, bloom2, bloom3, bloom4, cFP));

            break;
        }
    }
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomAlgorithm <KSIZE_1>;
template class DebloomAlgorithm <KSIZE_2>;
template class DebloomAlgorithm <KSIZE_3>;
template class DebloomAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
