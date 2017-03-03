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


#ifndef _GATB_CORE_BGLUE_ALGO_HPP_
#define _GATB_CORE_BGLUE_ALGO_HPP_

#include "unionFind.hpp"
#include <atomic>
#include <set>
#include <vector>
#include <string>
#include <mutex>
#include <unordered_map>
#include <BooPHF/BooPHF.h>
#include <ctime> // for time
#include <iostream> // for time (and maybe other things?)
#include <iomanip> // for cout mods
#include "ThreadPool.h"

#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/storage/impl/StorageTools.hpp>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/kmer/impl/PartiInfo.hpp>   // for repartitor 
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/collections/impl/BooPHF.hpp>


//heh at this point I could have maybe just included gatb_core.hpp but well, no circular dependencies, this file is part of gatb-core now.

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;



namespace gatb { namespace core { namespace debruijn { namespace impl  {

    template<size_t SPAN>
void bglue(Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int minSize, 
        int nb_threads, 
        int minimizer_type, 
        bool verbose
        );

}}}}

#endif

