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
        tools::collections::Iterable<T>* solidIterable,
        size_t                      kmerSize,
        tools::collections::impl::BloomFactory::Kind   bloomKind = tools::collections::impl::BloomFactory::CacheCoherent,
        const std::string&          debloomUri = "debloom",
        tools::misc::IProperties*   options    = 0
    );

    /** */
    ~DebloomAlgorithm ();

    /** */
    void execute ();

    /** Get the iterable over the computed critical FP kmers.
     * \return the cFP  kmers iterable. */
    tools::collections::Iterable<T>* getCriticalKmers ()  { return _criticalIterable; }

private:

    /** */
    virtual gatb::core::tools::collections::impl::Bloom<T>* createBloom ();

    void end_debloom_partition (gatb::core::tools::collections::impl::Hash16<T>& set, std::string& inputUri, std::string& outputUri);

    size_t       _kmerSize;
    tools::collections::impl::BloomFactory::Kind _bloomKind;
    std::string  _debloomUri;

    tools::collections::Iterable<T>* _solidIterable;
    void setSolidIterable (tools::collections::Iterable<T>* solidIterable)  {  _solidIterable = solidIterable; }

    /** */
    tools::collections::Iterable<T>* _criticalIterable;
    void setCriticalIterable (tools::collections::Iterable<T>* criticalIterable)
    {
        if (_criticalIterable != 0)  { delete _criticalIterable; }
        _criticalIterable = criticalIterable;
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _DEBLOOM_HPP_ */

