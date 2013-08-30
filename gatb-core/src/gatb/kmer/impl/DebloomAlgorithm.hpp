/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Debloom, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifndef _DEBLOOM_HPP_
#define _DEBLOOM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
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

template<typename T> class DebloomAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** */
    DebloomAlgorithm (
        tools::collections::impl::Product<tools::collections::impl::CollectionFile>& product,
        tools::collections::Iterable<T>* solidIterable,
        size_t                      kmerSize,
        tools::collections::impl::BloomFactory::Kind   bloomKind = tools::collections::impl::BloomFactory::CacheCoherent,
        const std::string&          debloomUri = "debloom",
        size_t                      max_memory = 1000,
        tools::misc::IProperties*   options    = 0
    );

    /** */
    ~DebloomAlgorithm ();

    /** */
    void execute ();

    /** Get the iterable over the computed critical FP kmers.
     * \return the cFP  kmers iterable. */
    tools::collections::Iterable<T>* getCriticalKmers ()  { return _criticalCollection->iterable(); }

private:

    /** */
    virtual gatb::core::tools::collections::impl::Bloom<T>* createBloom (tools::collections::Iterable<T>* solidIterable);

    void end_debloom_partition (
        gatb::core::tools::collections::impl::Hash16<T>& set,
        gatb::core::tools::dp::Iterator<T>*              inputIterator,
        gatb::core::tools::collections::Bag<T>*          outputBag
    );

    /** */
    tools::collections::impl::Product<tools::collections::impl::CollectionFile>& _product;

    size_t       _kmerSize;
    tools::collections::impl::BloomFactory::Kind _bloomKind;
    std::string  _debloomUri;
    size_t       _max_memory;

    tools::collections::Iterable<T>* _solidIterable;
    void setSolidIterable (tools::collections::Iterable<T>* solidIterable)  {  SP_SETATTR(solidIterable); }

    /** */
    tools::collections::impl::CollectionNode<tools::collections::impl::CollectionFile,T>* _criticalCollection;
    void setCriticalCollection (tools::collections::impl::CollectionNode<tools::collections::impl::CollectionFile,T>* criticalCollection)
    { _criticalCollection = criticalCollection; }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _DEBLOOM_HPP_ */

