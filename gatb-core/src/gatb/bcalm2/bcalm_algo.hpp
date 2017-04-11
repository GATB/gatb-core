/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014-2016  INRIA
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


#ifndef _GATB_CORE_BCALM_ALGO_HPP_
#define _GATB_CORE_BCALM_ALGO_HPP_

#include <gatb/tools/storage/impl/Storage.hpp>


namespace gatb { namespace core { namespace debruijn { namespace impl  {

    template<size_t SPAN>
void bcalm2(gatb::core::tools::storage::impl::Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int abundance, 
        int minSize, 
        int nb_threads, 
        int minimizer_type, 
        bool verbose
        );

}}}}

#endif
