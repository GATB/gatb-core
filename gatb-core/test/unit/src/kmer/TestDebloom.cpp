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
#include <gatb/kmer/impl/BloomAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>

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
class TestDebloom : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDebloom);

        CPPUNIT_TEST_GATB (Debloom_check1);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    () {}
    void tearDown () {}

        // SMALL VALUE NEEDED because continuous integration servers are not very powerful...
        static const u_int64_t MAX_MEMORY = 1000;

    /********************************************************************************/
    void Debloom_check1 ()
    {
        size_t kmerSize = 11;
        size_t miniSize = 8;
        size_t nks      = 1;

        const char* seqs[] = {
            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
        } ;

        /** We configure parameters for a SortingCountAlgorithm object. */
        IProperties* params = SortingCountAlgorithm<>::getDefaultProperties();  LOCAL (params);
        params->setInt (STR_KMER_SIZE,          kmerSize);
        params->setInt (STR_MINIMIZER_SIZE,     miniSize);
        params->setInt (STR_MAX_MEMORY,         MAX_MEMORY);
        params->setInt (STR_KMER_ABUNDANCE_MIN, nks);
        params->setStr (STR_URI_OUTPUT,         "foo");

        /** We create a SortingCountAlgorithm object. */
        SortingCountAlgorithm<> sortingCount (new BankStrings (seqs, ARRAY_SIZE(seqs)), params);

        /** We launch DSK. */
        sortingCount.execute();

        CPPUNIT_ASSERT (sortingCount.getSolidCounts()->getNbItems() == (int64_t)(strlen(seqs[0]) - kmerSize + 1) );

        /** We get the storage instance. */
        Storage* storage = sortingCount.getStorage();
        LOCAL (storage);

        Partition<SortingCountAlgorithm<>::Count>& counts = storage->getGroup("dsk").getPartition<Kmer<>::Count> ("solid");

        /** We create a bloom instance. */
        float nbitsPerKmer = DebloomAlgorithm<>::getNbBitsPerKmer (kmerSize, DEBLOOM_ORIGINAL);
        BloomAlgorithm<> bloom (*storage, &counts, kmerSize, nbitsPerKmer, 0, BLOOM_BASIC);
        bloom.execute ();

        /** We create a debloom instance. */
        DebloomAlgorithm<> debloom (
            storage->getGroup("bloom"),
            storage->getGroup("debloom"),
            &counts, kmerSize, miniSize, 1000, 0, BLOOM_BASIC, DEBLOOM_ORIGINAL
        );

        /** We launch the debloom. */
        debloom.execute();

        /** The following values have been computed with the original minia. */
        u_int64_t values[] =
        {
            0xc0620,    0x288f40,   0x188f40,   0x2aaa29,   0x8000b,    0x200881,   0x288081,   0x820db,    0x52e23,    0x2888f,
            0xaaa8b,    0x28838d,   0x20000,    0xa93ab,    0x2c18d,    0x2ba89,    0x183600,   0xea00b,    0x1a4ea0,   0xf8585
        };
        set<Kmer<>::Type> okValues;
        
        for (unsigned int i = 0; i < ARRAY_SIZE(values); i++)
        {
            Kmer<>::Type val; val.setVal(values[i]);
            okValues.insert(val);
        }

        CPPUNIT_ASSERT (debloom.getCriticalKmers()->getNbItems() == ARRAY_SIZE(values));

        /** We iterate the cFP kmers. */
        set<Kmer<>::Type> checkValues;
        Iterator<Kmer<>::Type>* iter = debloom.getCriticalKmers()->iterator();
        LOCAL (iter);

        for (iter->first(); !iter->isDone(); iter->next())
        {
            set<Kmer<>::Type>::iterator lookup = okValues.find (iter->item());
            CPPUNIT_ASSERT (lookup != okValues.end());

            checkValues.insert (iter->item());
        }

        CPPUNIT_ASSERT (checkValues.size() == okValues.size());
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestDebloom);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestDebloom);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

