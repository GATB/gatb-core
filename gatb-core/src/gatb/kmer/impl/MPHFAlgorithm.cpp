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

#include <gatb/kmer/impl/MPHFAlgorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <limits>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

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

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const char* messages[] = {
    "MPHF: initialization                   ",
    "MPHF: build hash function              ",
    "MPHF: assign values                    ",
    "MPHF: populate                         "
};

/** First tried to set the constant in the hpp file but got the following error:
 *  "error: a function call cannot appear in a constant-expression"
 *  Solved by putting it in the cpp...
 *      => http://stackoverflow.com/questions/2738435/using-numeric-limitsmax-in-constant-expressions
 */
template<size_t span, typename Abundance_t>
const Abundance_t MPHFAlgorithm<span,Abundance_t>::MAX_ABUNDANCE = std::numeric_limits<Abundance_t>::max();

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Abundance_t>
MPHFAlgorithm<span,Abundance_t>::MPHFAlgorithm (
    Group&              group,
    const std::string&  name,
    Iterable<Count>*    solidCounts,
    Iterable<Type>*     solidKmers,
    bool                buildOrLoad,
    IProperties*        options
)
    :  Algorithm("emphf", 1, options), _group(group), _name(name), _buildOrLoad(buildOrLoad),
       _solidCounts(0), _solidKmers(0), _map(0), _dataSize(0), _nb_abundances_above_precision(0), _progress(0)
{
    /** We keep a reference on the solid kmers. */
    setSolidCounts (solidCounts);

    /** We keep a reference on the solid kmers. */
    setSolidKmers (solidKmers);

    /** We build the hash object. */
    setMap (new Map());

    /** We gather some statistics. */
    getInfo()->add (1, "enabled", "%d", Map::enabled);

    /** In case of load, we load the mphf and populate right now. */
    if (Map::enabled == true && buildOrLoad == false)
    {
        /** We load the hash object from the dedicated storage group. */
        {   TIME_INFO (getTimeInfo(), "load");
            _map->load (_group, _name);
        }

        /** We populate the hash table. */
        populate ();
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span,typename Abundance_t>
MPHFAlgorithm<span,Abundance_t>::~MPHFAlgorithm ()
{
    /** Cleanup */
    setSolidCounts (0);
    setSolidKmers  (0);
    setMap         (0);
    setProgress    (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Abundance_t>
void MPHFAlgorithm<span,Abundance_t>::execute ()
{
    /** We check whether we can use such a type. */
    if (Map::enabled == true && _buildOrLoad == true)
    {
        size_t dataSize=0;

        /** We need a progress object. */
        tools::dp::IteratorListener* delegate = createIteratorListener(0,"");  LOCAL (delegate);
        setProgress (new ProgressCustom(delegate));

        /** We build the hash. */
        {   TIME_INFO (getTimeInfo(), "build");
            _map->build (*_solidKmers, _progress);
        }

        /** We save the hash object in the dedicated storage group. */
        {   TIME_INFO (getTimeInfo(), "save");
            _dataSize = _map->save (_group, _name);
        }

        /** We populate the hash table. */
        populate ();
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span,typename Abundance_t>
float MPHFAlgorithm<span,Abundance_t>::getNbBitsPerKmer () const
{
    float nbitsPerKmer = sizeof(Abundance_t)*8;
    return nbitsPerKmer;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span,typename Abundance_t>
void MPHFAlgorithm<span,Abundance_t>::populate ()
{
    size_t nb_iterated = 0;
    size_t n = _map->size();

    _nb_abundances_above_precision = 0;

    /** We need a progress object. */
    tools::dp::IteratorListener* delegate = createIteratorListener(_solidCounts->getNbItems(),messages[3]);  LOCAL (delegate);
    setProgress (new ProgressCustom(delegate));

    SubjectIterator<Count>* itKmers = new SubjectIterator<Count> (_solidCounts->iterator(), _solidCounts->getNbItems()/100);
    itKmers->addObserver (_progress);
    LOCAL (itKmers);

    // set counts and at the same time, test the mphf
    for (itKmers->first(); !itKmers->isDone(); itKmers->next())
    {
        /** We get the hash code of the current item. */
        typename Map::Hash::Code h = _map->getCode (itKmers->item().value);

        /** Little check. */
        if (h >= n) {  throw Exception ("MPHF check: value out of bounds"); }

        /** We get the abundance of the current kmer. */
        int abundance = itKmers->item().abundance;

        if (abundance > MAX_ABUNDANCE)
        {
            _nb_abundances_above_precision++;
            abundance = MAX_ABUNDANCE;
        }

        /** We set the abundance of the current kmer. */
        _map->at (h) = abundance;

        nb_iterated ++;
    }

    if (nb_iterated != n && n > 3)
    {
        throw Exception ("ERROR during abundance population: itKmers iterated over %d/%d kmers only", nb_iterated, n);
    }

#if 1
    // you know what? let's test if the MPHF does not have collisions, it won't hurt.
    check ();
#endif

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "nb_keys",               "%ld",  _map->size());
    getInfo()->add (2, "data_size",             "%ld",  _dataSize);
    getInfo()->add (2, "bits_per_key",          "%.3f", (float)(_dataSize*8)/(float)_map->size());
    getInfo()->add (2, "prec",                  "%d",   MAX_ABUNDANCE);
    getInfo()->add (2, "nb_abund_above_prec",   "%d",   _nb_abundances_above_precision);
    getInfo()->add (1, getTimeInfo().getProperties("time"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span,typename Abundance_t>
void MPHFAlgorithm<span,Abundance_t>::check ()
{
    size_t nb_iterated = 0;

    Iterator<Count>* itKmers = _solidCounts->iterator();  LOCAL (itKmers);

    for (itKmers->first(); !itKmers->isDone(); itKmers->next())
    {
        Count& count = itKmers->item();

        /** We get the current abundance. */
        Abundance_t abundance = (*_map)[count.value];

        if (abundance!=count.abundance && abundance<MAX_ABUNDANCE)  {  throw Exception ("ERROR: MPHF isn't injective (abundance population failed)");  }

        nb_iterated ++;
    }

    if (nb_iterated != _map->size() && _map->size() > 3)
    {
        throw Exception ("ERROR during abundance population: itKmers iterated over %d/%d kmers only", nb_iterated, _map->size());
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span,typename Abundance_t>
MPHFAlgorithm<span,Abundance_t>::ProgressCustom::ProgressCustom (tools::dp::IteratorListener* ref)
  : tools::misc::impl::ProgressProxy (ref), nbReset(0)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span,typename Abundance_t>
void MPHFAlgorithm<span,Abundance_t>::ProgressCustom::reset (u_int64_t ntasks)
{
    const char* label = nbReset < sizeof(messages)/sizeof(messages[0]) ? messages[nbReset++] : "other";

    getRef()->reset(ntasks);
    getRef()->setMessage (label);
    getRef()->init();
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class MPHFAlgorithm <KSIZE_1>;
template class MPHFAlgorithm <KSIZE_2>;
template class MPHFAlgorithm <KSIZE_3>;
template class MPHFAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
