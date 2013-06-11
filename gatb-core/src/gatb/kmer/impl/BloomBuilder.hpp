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
#include <gatb/tools/misc/api/IProperty.hpp>

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
class BloomBuilder
{
public:

    /** */
    BloomBuilder (tools::dp::Iterator<kmer_type>* itKmers, size_t bloomSize, size_t nbHash, size_t nbCores=0);

    /** */
    ~BloomBuilder ();

    /** */
    tools::collections::impl::Bloom<kmer_type>*  build (tools::misc::IProperties* stats = 0);

private:

    tools::dp::Iterator<kmer_type>* _itKmers;
    void setItKmers (tools::dp::Iterator<kmer_type>* itKmers)  { SP_SETATTR(itKmers); }

    size_t _bloomSize;
    size_t _nbHash;
    size_t _nbCores;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_ */
