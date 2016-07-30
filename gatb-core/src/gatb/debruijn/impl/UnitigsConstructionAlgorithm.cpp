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

#include <gatb/debruijn/impl/UnitigsConstructionAlgorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/bcalm2/bcalm_algo.hpp>
#include <gatb/bcalm2/bglue_algo.hpp>

#include <queue>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace debruijn  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Graph: building unitigs using bcalm2  ";
static const char* progressFormat2 = "Graph: nb unitigs constructed : %-9d  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
UnitigsConstructionAlgorithm<span>::UnitigsConstructionAlgorithm (
    tools::storage::impl::Storage& storage,
    std::string                 unitigs_filename,
    size_t                      nb_cores,
    tools::misc::IProperties*   options
)
    : Algorithm("bcalm2-wrapper", nb_cores, options), _storage(storage), unitigs_filename(unitigs_filename)
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
template<size_t span>
UnitigsConstructionAlgorithm<span>::~UnitigsConstructionAlgorithm ()
{
}

/*********************************************************************/
/*********************************************************************/

template <size_t span>
void UnitigsConstructionAlgorithm<span>::execute ()
{
    int kmerSize =
            getInput()->getInt(STR_KMER_SIZE);
    int abundance = 
            getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
    int minimizerSize =
        getInput()->getInt(STR_MINIMIZER_SIZE);
    int nb_threads =
        getInput()->getInt(STR_NB_CORES);
    int minimizer_type =
        getInput()->getInt(STR_MINIMIZER_TYPE);
    bool verbose = getInput()->getInt(STR_VERBOSE);
        
    unsigned int nbThreads = this->getDispatcher()->getExecutionUnitsNumber();
    if ((unsigned int)nb_threads > nbThreads)
    {
        std::cout << "Uh. Unitigs graph construction called with nb_threads " << nb_threads << " but dispatcher has nbThreads " << nbThreads << std::endl;

    }

    bcalm2<span>(&_storage, unitigs_filename, kmerSize, abundance, minimizerSize, nbThreads, minimizer_type, verbose); 
    bglue<span> (&_storage, unitigs_filename, kmerSize,            minimizerSize, nbThreads, minimizer_type, verbose);

    /** We gather some statistics. */
    //getInfo()->add (1, "stats");
    //getInfo()->add (2, "nb_unitigs", "%ld", /* */);
    
    //getInfo()->add (1, "time");
    //getInfo()->add (2, "build", "%.3f", /* */);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
