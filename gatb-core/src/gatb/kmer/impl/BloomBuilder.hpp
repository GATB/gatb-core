/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file BloomBuilder.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Build bloom filter from an iterator of kmers
 */

#ifndef _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_
#define _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/collections/impl/Bloom.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/math/NativeInt8.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** \brief tbd
 *
 * This class is a builder of a bloom filter in which we insert kmers.
 */
template<size_t span> class BloomBuilder
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    /** */
    BloomBuilder (
        u_int64_t   bloomSize,
        size_t      nbHash,
        tools::collections::impl::BloomFactory::Kind   bloomKind = tools::collections::impl::BloomFactory::CacheCoherent,
        size_t      nbCores = 0
    )
        : _bloomSize (bloomSize), _nbHash (nbHash), _nbCores(nbCores), _bloomKind(bloomKind)
    {
    }

    /** */
    tools::collections::impl::Bloom<Type>*  build (
        tools::dp::Iterator<Count>* itKmers,
        tools::misc::IProperties* stats=0
    )
    {
        tools::misc::impl::TimeInfo ti;
        TIME_INFO (ti, "build_kmers_bloom");

        LOCAL (itKmers);

        /** We instantiate the bloom object. */
        tools::collections::impl::Bloom<Type>* bloom =
            tools::collections::impl::BloomFactory::singleton().createBloom<Type> (_bloomKind, _bloomSize, _nbHash);

        /** We launch the bloom fill. */
        tools::dp::impl::Dispatcher(_nbCores).iterate (itKmers,  BuildKmerBloom (*bloom));

        /** We gather some statistics. */
        if (stats != 0)
        {
            //stats->add (0, "bloom");
            stats->add (0, "filter_size", "%lld", _bloomSize);
            stats->add (0, "nb_hash_fct", "%d",   _nbHash);
        }

        /** We return the created bloom filter. */
        return bloom;
    }

    /** */
    tools::collections::impl::Bloom<Type>*  load (
        tools::collections::Iterable<tools::math::NativeInt8>* bloomIterable,
        tools::misc::IProperties* stats = 0)
    {
        tools::misc::impl::TimeInfo ti;
        TIME_INFO (ti, "load_bloom");

        LOCAL (bloomIterable);

        /** We instantiate the bloom object. */
        tools::collections::impl::Bloom<Type>* bloom =
            tools::collections::impl::BloomFactory::singleton().createBloom<Type> (_bloomKind, _bloomSize, _nbHash);

        /** We set the bloom with the provided array given as an iterable of NativeInt8 objects. */
        bloomIterable->getItems ((tools::math::NativeInt8*&)bloom->getArray());

        /** We gather some statistics. */
        if (stats != 0)
        {
            stats->add (0, "bloom");
            stats->add (1, "filter_size", "%lld", _bloomSize);
            stats->add (1, "nb_hash_fct", "%d",   _nbHash);
        }

        /** We return the created bloom filter. */
        return bloom;
    }


private:

    u_int64_t _bloomSize;
    size_t    _nbHash;
    size_t    _nbCores;
    tools::collections::impl::BloomFactory::Kind _bloomKind;

    /********************************************************************************/
    class BuildKmerBloom
    {
    public:
        void operator() (const Count& kmer)  {  _bloom.insert(kmer.value); }
        BuildKmerBloom (tools::collections::impl::Bloom<Type>& bloom)  : _bloom(bloom) {}
        tools::collections::impl::Bloom<Type>& _bloom;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_ */
