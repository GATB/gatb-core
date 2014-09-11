#ifdef WITH_MPHF
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
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/NativeInt8.hpp>

#include <gatb/tools/storage/impl/StorageTools.hpp>



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

using namespace emphf;

#define DEBUG(a)  printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "MPHF: read solid kmers              ";
static const char* progressFormat2 = "MPHF: build mphf                    ";
static const char* progressFormat3 = "MPHF: finalization                  ";

template<size_t span>
MPHFAlgorithm<span>::MPHFAlgorithm (
    Storage& storage,
    Iterable<Count>*    solidIterable,
    size_t              kmerSize,
    const std::string&  mphfUri,
    IProperties*        options,
    bool load_emphf
)
    :  Algorithm("emphf", 1, options), _storage(storage), _group(storage().getGroup ("mphf")),
       _kmerSize(kmerSize),
       _mphfUri(mphfUri),
       _solidIterable(0)
{
    /** We get a group for MPHF. */
    Group& group = _storage().getGroup ("mphf");

    setSolidIterable    (solidIterable);

    // TODO: should later be replaced by better serialization
    // and a constructor taking only _storage, like DebloomAlgorithm
    if (load_emphf)
    {

        Model model (kmerSize);
        long n = solidIterable->getNbItems();

        mphf_class = new MPHF(_kmerSize, n);

        mphf_class->setEMPHF(StorageTools::singleton().loadEMPHF<typename MPHF::BaseHasher>(_group, _mphfUri));

        mphf_class->populateAbundances(_solidIterable);

        //fprintf(stderr, "Loaded MPHF from disk, number of elements: %d\n", mphf_class->getSize());

    }
}

template<size_t span>
MPHFAlgorithm<span>::~MPHFAlgorithm ()
{
    setSolidIterable      (0);
}



template<size_t span>
void MPHFAlgorithm<span>::execute ()
{
    Model model (_kmerSize);
    long n = _solidIterable->getNbItems();

    mphf_class = new MPHF(_kmerSize, n, _solidIterable);

    mphf_class->populateAbundances(_solidIterable);

    // save the produced mphf
    StorageTools::singleton().saveEMPHF<typename MPHF::BaseHasher>(_group, _mphfUri, static_cast< void* >( &( mphf_class->mphf) ));
}

template<size_t span>
float MPHFAlgorithm<span>::getNbBitsPerKmer () const
{
    float nbitsPerKmer = sizeof(mphf_abundance_t)*8;

    return nbitsPerKmer;
}

template<size_t span>
typename MPHFAlgorithm<span>::MPHF * MPHFAlgorithm<span>::getMPHF () 
{
    return mphf_class;
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
#endif
