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

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>

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

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

template <typename Item>
class BuildKmerExtension : public IteratorFunctor
{
public:
    void operator() (const Item& kmer)
    {
        /** We configure the neighbor kmers iterator for a given source kmer. */
        _itNeighbors.setSource (kmer);

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

    ~BuildKmerExtension ()  {  _extendBag.flush();  }

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
template<typename T>
DebloomAlgorithm<T>::DebloomAlgorithm (
    tools::collections::Iterable<T>* solidIterable,
    size_t              kmerSize,
    BloomFactory::Kind  bloomKind,
    const std::string&  debloomUri,
    IProperties*        options
)
    :  Algorithm("debloom", options), _kmerSize(kmerSize), _bloomKind(bloomKind), _debloomUri("debloom"),
       _solidIterable(0), _criticalIterable(0)
{
    setSolidIterable    (solidIterable);
    setCriticalIterable (new IterableFile<T> (_debloomUri, 2000));
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
DebloomAlgorithm<T>::~DebloomAlgorithm ()
{
    setSolidIterable    (0);
    setCriticalIterable (0);
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
void DebloomAlgorithm<T>::execute ()
{
    Model<T> model (_kmerSize);

    /*************************************************/
    /** We create an iterator over the solid kmers.  */
    /*************************************************/
#if 1
    Iterator<T>* itKmers = createIterator<T> (
        _solidIterable->iterator(),
        _solidIterable->getNbItems(),
        "iterate solid kmers"
    );
#else
    Iterator<T>* itKmers = _solidIterable->iterator();
#endif
    LOCAL (itKmers);

    /*************************************************/
    /** We create the debloom file.                 */
    /*************************************************/

    /** First, we delete the debloom file if already existing. */
    System::file().remove (_debloomUri);

    Bag<T>* debloomFile = new BagFile<T>(_debloomUri);
    LOCAL (debloomFile);

    /*************************************************/
    /** We fill the debloom file.                    */
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "fill debloom file");

        /** We create a bloom with inserted solid kmers. */
        Bloom<T>* bloom = createBloom ();
        LOCAL (bloom);

        /** We iterate the kmers and add them into the bloom filter. */
        getDispatcher()->iterate (itKmers, BuildKmerExtension<T> (model, bloom, debloomFile));
    }

    /** We make sure everything is put into the extension file. */
    debloomFile->flush();

    /*************************************************/
    /** We compute the final cFP file.               */
    /*************************************************/

    size_t max_memory = 1000; //getInput()->getInt (DSK::STR_MAX_MEMORY);

    Hash16<T> partition (max_memory);

    string inputUri  = _debloomUri;
    string outputUri = _debloomUri + "2";

    {
        TIME_INFO (getTimeInfo(), "finalize debloom file");

        /** We iterate the solid kmers. */
        for (itKmers->first(); !itKmers->isDone(); itKmers->next())
        {
            /** We add the current kmer into the partition. */
            partition.insert (itKmers->item());

            /** We may have reach the maximum number of items in the partition. */
            if (partition.size() >= partition.getMaxNbItems())
            {
                /** We exclude the partition content from the critical false positive file. */
                end_debloom_partition (partition, inputUri, outputUri);
            }
        }

        /** We exclude the partition content from the critical false positive file. */
        end_debloom_partition (partition, inputUri, outputUri);
    }

    /** We swap the filenames (because the last 'end_debloom_partition' call made one too much). */
    std::swap (inputUri, outputUri);

    /** We make sure that 1) the final cFP file has the good name  2) we remove the other temporary debloom file. */
    if (outputUri != _debloomUri)
    {
        /** We rename the temporary file to the final name. */
        System::file().rename (outputUri, inputUri);
    }
    else
    {
        /** We delete the temporary file. */
        System::file().remove (inputUri);
    }
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

    ~EndDebloomPartition ()
    {
        /** We need to flush the cache to be sure to get all the items in the actual cached bag. */
        _output.flush();
    }

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
template<typename T>
void DebloomAlgorithm<T>::end_debloom_partition (Hash16<T>& partition, string& inputUri, string& outputUri)
{
    /** First, we delete the output debloom file if already existing. */
    System::file().remove (outputUri);

    /** We need an input and an output. */
    IteratorFile<T>* debloomInput  = new IteratorFile <T> (inputUri);
    Bag<T>*          debloomOutput = new BagFile      <T> (outputUri);

    LOCAL (debloomInput);
    LOCAL (debloomOutput);

    /** We dispatch the iteration of the current cFP bag. */
    getDispatcher()->iterate (debloomInput, EndDebloomPartition<T> (partition, debloomOutput));

    /** We make sure we output all the items. */
    debloomOutput->flush ();

    /** We swap the filenames. */
    std::swap (inputUri, outputUri);

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
template<typename T>
Bloom<T>* DebloomAlgorithm<T>::createBloom ()
{
    TIME_INFO (getTimeInfo(), "fill bloom filter");

    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

    u_int64_t solidKmersNb = _solidIterable->getNbItems();

    u_int64_t estimatedBloomSize = (u_int64_t) (solidKmersNb * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator<T>* itKmers = createIterator<T> (
        _solidIterable->iterator(),
        solidKmersNb,
        "fill bloom filter  "
    );
    LOCAL (itKmers);

    /** We use a bloom builder. */
    BloomBuilder<T> builder (itKmers, estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER), _bloomKind, getDispatcher()->getExecutionUnitsNumber());

    /** We instantiate the bloom object. */
    IProperties* bloomProps = new Properties();
    Bloom<T>* bloom = builder.build (bloomProps);

    getInfo()->add (1, bloomProps);
    getInfo()->add (2, "nbits per kmer", "%f", NBITS_PER_KMER);

    /** We return the created bloom filter. */
    return bloom;
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomAlgorithm <gatb::core::tools::math::NativeInt64>;
#ifdef INT128_FOUND
template class DebloomAlgorithm <gatb::core::tools::math::NativeInt128>;
#else
template class DebloomAlgorithm <gatb::core::tools::math::LargeInt<2> >;
#endif

template class DebloomAlgorithm <gatb::core::tools::math::LargeInt<3> >;
template class DebloomAlgorithm <gatb::core::tools::math::LargeInt<4> >;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
