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

/** \file BloomBuilder.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Build bloom filter from an iterator of kmers
 */

#ifndef _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_
#define _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/collections/api/Iterable.hpp>
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
template<size_t span=KMER_DEFAULT_SPAN> class BloomBuilder
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    /** */
    BloomBuilder (
        u_int64_t   bloomSize,
        size_t      nbHash,
        tools::misc::BloomKind bloomKind = tools::misc::BLOOM_DEFAULT,
        size_t      nbCores = 0,
		int		  min_abundance =0
    )
        : _bloomSize (bloomSize), _nbHash (nbHash), _nbCores(nbCores), _bloomKind(bloomKind), _min_abundance(min_abundance)
    {
    }

    /** */
    tools::collections::impl::IBloom<Type>*  build (
        tools::dp::Iterator<Count>* itKmers,
        tools::misc::IProperties* stats=0
    )
    {
        tools::misc::impl::TimeInfo ti;
        TIME_INFO (ti, "build_kmers_bloom");

        LOCAL (itKmers);

        /** We instantiate the bloom object. */
        tools::collections::impl::IBloom<Type>* bloom =
            tools::collections::impl::BloomFactory::singleton().createBloom<Type> (_bloomKind, _bloomSize, _nbHash);

        /** We launch the bloom fill. */
        tools::dp::impl::Dispatcher(_nbCores).iterate (itKmers,  BuildKmerBloom (*bloom,_min_abundance));

        /** We gather some statistics. */
        if (stats != 0)
        {
            //stats->add (0, "bloom");
            stats->add (0, "size",    "%lld", _bloomSize);
            stats->add (0, "nb_hash", "%d",   _nbHash);
        }

        /** We return the created bloom filter. */
        return bloom;
    }

    /** */
    tools::collections::impl::IBloom<Type>*  load (
        tools::collections::Iterable<tools::math::NativeInt8>* bloomIterable,
        tools::misc::IProperties* stats = 0)
    {
        tools::misc::impl::TimeInfo ti;
        TIME_INFO (ti, "load_bloom");

        LOCAL (bloomIterable);

        /** We instantiate the bloom object. */
        tools::collections::impl::IBloom<Type>* bloom =
            tools::collections::impl::BloomFactory::singleton().createBloom<Type> (_bloomKind, _bloomSize, _nbHash);

        /** We set the bloom with the provided array given as an iterable of NativeInt8 objects. */
        bloomIterable->getItems ((tools::math::NativeInt8*&)bloom->getArray());

        /** We gather some statistics. */
        if (stats != 0)
        {
            stats->add (0, "bloom", "");
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
	int _min_abundance;
	
    tools::misc::BloomKind _bloomKind;

    /********************************************************************************/
    class BuildKmerBloom
    {
    public:
        void operator() (const Count& kmer)  { if(kmer.abundance >= _min_abundance)  _bloom.insert(kmer.value); }
        BuildKmerBloom (tools::collections::impl::IBloom<Type>& bloom, int min_abundance=0)  : _bloom(bloom),_min_abundance(min_abundance)  {}
        tools::collections::impl::IBloom<Type>& _bloom;
		int _min_abundance;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_ */
