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


#ifndef _GATB_CORE_LINK_TIGS_HPP_
#define _GATB_CORE_LINK_TIGS_HPP_

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bank/api/IBank.hpp>

namespace gatb { namespace core { namespace debruijn { namespace impl  {


    template<size_t SPAN>
    void link_tigs( std::string prefix, int kmerSize, int nb_threads, uint64_t &nb_unitigs, bool verbose);

    template<size_t span>
    void link_unitigs_pass(const std::string unitigs_filename, bool verbose, const int pass, const int kmerSize);
    
}}}}

#endif
