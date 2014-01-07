/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Debloom, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file DebloomAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Debloom algorithm, ie. compute false positive sets for a Bloom filter
 */

#ifndef _DEBLOOM_HPP_
#define _DEBLOOM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/Hash16.hpp>

#include <gatb/tools/collections/impl/Product.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<typename ProductFactory, size_t span=KMER_DEFAULT_SPAN>
class DebloomAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::Model Model;
    typedef typename kmer::impl::Kmer<span>::Type  Type;
    typedef typename kmer::impl::Kmer<span>::Count Count;

    /** */
    DebloomAlgorithm (
        tools::collections::impl::Product<ProductFactory>& product,
        tools::collections::Iterable<Count>* solidIterable,
        size_t                      kmerSize,
        size_t                      max_memory = 1000,
        size_t                      nb_cores   = 0,
        tools::collections::impl::BloomFactory::Kind   bloomKind = tools::collections::impl::BloomFactory::CacheCoherent,
        const std::string&          debloomUri = "debloom",
        tools::misc::IProperties*   options    = 0
    );

    /** */
    ~DebloomAlgorithm ();

    /** */
    void execute ();

    /** Get the collection for the computed critical FP kmers.
     * \return the cFP  kmers collection. */
    tools::collections::Collection<Type>* getCriticalKmers ()  { return _criticalCollection; }

private:

    /** */
    virtual gatb::core::tools::collections::impl::Bloom<Type>* createBloom (
        tools::collections::Iterable<Count>* solidIterable,
        tools::misc::IProperties* props
    );

    void end_debloom_partition (
        gatb::core::tools::collections::impl::Hash16<Type>& set,
        gatb::core::tools::dp::Iterator<Type>*              inputIterator,
        gatb::core::tools::collections::Bag<Type>*          outputBag
    );

    /** */
    tools::collections::impl::Product<ProductFactory>& _product;

    size_t       _kmerSize;
    tools::collections::impl::BloomFactory::Kind _bloomKind;
    std::string  _debloomUri;
    size_t       _max_memory;

    tools::collections::Iterable<Count>* _solidIterable;
    void setSolidIterable (tools::collections::Iterable<Count>* solidIterable)  {  SP_SETATTR(solidIterable); }

    /** */
    tools::collections::impl::CollectionNode<Type>* _criticalCollection;
    void setCriticalCollection (tools::collections::impl::CollectionNode<Type>* criticalCollection)
    { _criticalCollection = criticalCollection; }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _DEBLOOM_HPP_ */

