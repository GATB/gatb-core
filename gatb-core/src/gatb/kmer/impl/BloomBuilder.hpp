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
template<typename T> class BloomBuilder
{
public:

    /** */
    BloomBuilder (tools::dp::Iterator<T>* itKmers, size_t bloomSize, size_t nbHash, size_t nbCores=0)
        : _itKmers (0), _bloomSize (bloomSize), _nbHash (nbHash), _nbCores(nbCores)
    {
        setItKmers (itKmers);
    }

    /** */
    ~BloomBuilder ()  {  setItKmers (0);  }

    /** */
    tools::collections::impl::Bloom<T>*  build (tools::misc::IProperties* stats = 0)
    {
        tools::misc::impl::TimeInfo ti;
        TIME_INFO (ti, "build kmers bloom");

        /** We instantiate the bloom object. */
        tools::collections::impl::Bloom<T>* bloom = new tools::collections::impl::BloomCacheCoherent<T> (_bloomSize, _nbHash);
      //  tools::collections::impl::Bloom<T>* bloom = new tools::collections::impl::BloomSynchronized<T> (_bloomSize, _nbHash);

        /** We launch the bloom fill. */
        tools::dp::impl::ParallelCommandDispatcher(_nbCores).iterate (_itKmers,  BuildKmerBloom (*bloom));

        /** We gather some statistics. */
        if (stats != 0)
        {
            stats->add (0, "bloom");
            stats->add (1, "filter size", "%d", _bloomSize);
            stats->add (1, "nb hash fct", "%d", _nbHash);
        }

        /** We return the created bloom filter. */
        return bloom;
    }

private:

    tools::dp::Iterator<T>* _itKmers;
    void setItKmers (tools::dp::Iterator<T>* itKmers)  { SP_SETATTR(itKmers); }

    size_t _bloomSize;
    size_t _nbHash;
    size_t _nbCores;

    /********************************************************************************/
    class BuildKmerBloom : public tools::dp::impl::IteratorFunctor
    {
    public:
        void operator() (const T& kmer)  {  _bloom.insert(kmer); }
        BuildKmerBloom (tools::collections::impl::Bloom<T>& bloom)  : _bloom(bloom) {}
        tools::collections::impl::Bloom<T>& _bloom;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_ */
