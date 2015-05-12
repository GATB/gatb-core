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

/** \file DebloomAlgorithmFactory.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Factory for Debloom algorithms
 */

#ifndef _DEBLOOM_ALGORITHM_FACTORY_HPP_
#define _DEBLOOM_ALGORITHM_FACTORY_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/DebloomAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomMinimizerAlgorithm.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class DebloomAlgorithmFactory
{
public:

    typedef typename kmer::impl::Kmer<span>::Count Count;

    static DebloomAlgorithm<span>* create (
        tools::misc::DebloomImpl impl,
        tools::storage::impl::Group&    bloomGroup,
        tools::storage::impl::Group&    debloomGroup,
        tools::storage::impl::Partition<Count>* solidIterable,
        size_t                      kmerSize,
        size_t                      miniSize,
        size_t                      max_memory = 0,
        size_t                      nb_cores   = 0,
        tools::misc::BloomKind      bloomKind     = tools::misc::BLOOM_DEFAULT,
        tools::misc::DebloomKind    cascadingKind = tools::misc::DEBLOOM_DEFAULT,
        const std::string&          debloomUri = "debloom",
        tools::misc::IProperties*   options    = 0,
        tools::storage::impl::Group*    minimizersGroup = 0
    )
    {
        switch (impl)
        {
            case tools::misc::DEBLOOM_IMPL_BASIC:
                return new DebloomAlgorithm<span> (
                    bloomGroup, debloomGroup, solidIterable, kmerSize, miniSize, max_memory, nb_cores,
                    bloomKind, cascadingKind, debloomUri, options
                );

            case tools::misc::DEBLOOM_IMPL_MINIMIZER:
            case tools::misc::DEBLOOM_IMPL_DEFAULT:
            default:
                return new DebloomMinimizerAlgorithm<span> (
                    bloomGroup, debloomGroup, solidIterable, kmerSize, miniSize, max_memory, nb_cores,
                    bloomKind, cascadingKind, debloomUri, options, minimizersGroup
                );
        }
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _DEBLOOM_ALGORITHM_FACTORY_HPP_ */

