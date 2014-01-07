/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file gatb_core.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Global header file
 */

#ifndef _GATB_CORE_HPP_
#define _GATB_CORE_HPP_

/********************************************************************************/

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/ContainerSet.hpp>
#include <gatb/tools/collections/impl/Hash16.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/collections/impl/Product.hpp>
#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankStrings.hpp>
#include <gatb/bank/impl/BankSplitter.hpp>
#include <gatb/bank/impl/BankRandom.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/GraphHelpers.hpp>

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;

/********************************************************************************/

#endif /* _GATB_CORE_HPP_ */
