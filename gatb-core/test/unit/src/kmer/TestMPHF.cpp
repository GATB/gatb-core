/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

#include <CppunitCommon.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/BankStrings.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/MPHFAlgorithm.hpp>

#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestMPHF : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestMPHF);

        CPPUNIT_TEST_GATB (MPHF_check1);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    () {}
    void tearDown () {}


    /********************************************************************************/
    void MPHF_check1 ()
    {
        size_t kmerSize = 11;
        size_t nks      = 1;

        const char* seqs[] = {
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
        } ;

        /** We create a storage instance. */
        Storage storage (STORAGE_FILE, "test");

        /** We create a DSK instance. */
        SortingCountAlgorithm<> sortingCount (&storage, new BankStrings (seqs, ARRAY_SIZE(seqs)), kmerSize, nks);

        /** We launch DSK. */
        sortingCount.execute();

        CPPUNIT_ASSERT (sortingCount.getSolidKmers()->getNbItems() == (strlen(seqs[0]) - kmerSize + 1) );

        /** We create a mphf instance. */
        MPHFAlgorithm<> mphf (storage, sortingCount.getSolidKmers(), kmerSize);

        /** We actually execute the mphf construction. */
        mphf.execute();

        MPHFAlgorithm<>::MPHF *mphf_class = mphf.getMPHF();

        // below are quick tests
        
        CPPUNIT_ASSERT (mphf_class->getSize() == 130);
       
        typedef typename Kmer<32>::Model Model;
        typedef typename Kmer<32>::Type  Type;
        Model model (11);

        Type kmer = model.codeSeed ("ACCATGTATAA", Data::ASCII);
        mphf_class->set(kmer,4);

        unsigned char val = 0;
        mphf_class->get(kmer,&val);

        CPPUNIT_ASSERT (val == 4);
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestMPHF);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestMPHF);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

