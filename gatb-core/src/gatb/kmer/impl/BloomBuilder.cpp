/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <math.h>

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/
namespace gatb  { namespace core  {  namespace kmer  {  namespace impl {
/********************************************************************************/

/********************************************************************************/
struct BuildKmerBloom
{
    void operator() (const kmer_type& kmer)  {  _bloom.insert(kmer); }
    BuildKmerBloom (Bloom<kmer_type>& bloom)  : _bloom(bloom) {}
    Bloom<kmer_type>& _bloom;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BloomBuilder::BloomBuilder (Iterator<kmer_type>* itKmers, size_t bloomSize, size_t nbHash, size_t nbCores)
    : _itKmers (0), _bloomSize (bloomSize), _nbHash (nbHash), _nbCores(nbCores)
{
    setItKmers (itKmers);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BloomBuilder::~BloomBuilder ()
{
    setItKmers (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bloom<kmer_type>* BloomBuilder::build (IProperties* stats)
{
    TimeInfo ti;
    TIME_INFO (ti, "build kmers bloom");

    /** We instantiate the bloom object. */
    Bloom<kmer_type>* bloom = new BloomSynchronized<kmer_type> (_bloomSize, _nbHash);

    /** We launch the bloom fill. */
    ParallelCommandDispatcher(_nbCores).iterate (*_itKmers,  BuildKmerBloom (*bloom));

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

/********************************************************************************/
} } } }
/********************************************************************************/
