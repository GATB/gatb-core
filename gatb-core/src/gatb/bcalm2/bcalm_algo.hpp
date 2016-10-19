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

#include <assert.h>
#include <iostream>
#include <memory>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <tuple>
#include "ograph.h"

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

#include <thread>
#include <atomic>
#include "lockstdqueue.h"

#include "ThreadPool.h"


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




#define get_wtime() chrono::system_clock::now()
#ifndef diff_wtime
#define diff_wtime(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count()
#endif

//#define BINSEQ // "graph4 is not ready" according to antoine. also, initBinSeq provokes segfault at end of bcalm

#ifdef BINSEQ
#include "binSeq.h"
#define BUCKET_STR_TYPE binSeq
#define TO_BUCKET_STR(x) binSeq(x)
#define FROM_BUCKET_STR(x) (x.str())
#else
#define BUCKET_STR_TYPE string
#define TO_BUCKET_STR(x) x
#define FROM_BUCKET_STR(x) x
#endif


// timing-related variables

#define THREAD_SAFE_TIMING
#ifdef THREAD_SAFE_TIMING
typedef std::atomic<double> atomic_double;
#else
#define atomic_double_add(d1,d2) d1 += d2;
typedef double atomic_double;
#endif


namespace gatb { namespace core { namespace debruijn { namespace impl  {

    template<size_t SPAN>
void bcalm2(Storage* storage, 
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
