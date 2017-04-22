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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_UTIGS_CONSTR_ALGORITHM_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_UTIGS_CONSTR_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/kmer/impl/Model.hpp> // for KMER_DEFAULT_SPAN and so on
                                                                                                                                                                                
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/bcalm2/bcalm_algo.hpp>


/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

/** \brief Computation of the unitigs of a Graph, using BCALM 2
 *
 */
template <size_t span=KMER_DEFAULT_SPAN>
class UnitigsConstructionAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:
    /** Constructor.
     * \param[in] graph : graph from which we look for branching nodes
     * \param[in] nb_cores : number of cores to be used; 0 means all available cores
     * \param[in] options : extra options
     */
    UnitigsConstructionAlgorithm (
        tools::storage::impl::Storage& storage,
        std::string                 unitigs_filename,
        size_t                      nb_cores = 0,
        tools::misc::IProperties*   options  = 0,
        bool do_bcalm = true,
        bool do_bglue = true,
        bool do_links = true
    );
    
    /** Destructor. */
    ~UnitigsConstructionAlgorithm ();

    /** Get an option parser for branching parameters. Dynamic allocation, so must be released when no more used.
     * \return an instance of IOptionsParser.
     */
    static tools::misc::IOptionsParser* getOptionsParser ();

    /** \copydoc tools::misc::impl::Algorithm::execute */
    void execute ();
    
    int kmerSize;

    uint64_t nb_unitigs;

private:

    tools::storage::impl::Storage& _storage;
    std::string unitigs_filename;
    bool do_bcalm, do_bglue, do_links;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif

