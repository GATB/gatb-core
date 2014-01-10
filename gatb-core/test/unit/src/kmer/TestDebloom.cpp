/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
 *****************************************************************************/

#include <CppunitCommon.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/BankStrings.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
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


    /********************************************************************************/
    void Debloom_check1 ()
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

        /** We create a debloom instance. */
        DebloomAlgorithm<> debloom (storage, sortingCount.getSolidKmers(), kmerSize, 1000, 0, BloomFactory::Synchronized);

        /** We launch the debloom. */
        debloom.execute();

        /** The following values have been computed with the original minia. */
        u_int64_t values[] =
        {
            0xc0620,    0x288f40,   0x188f40,   0x2aaa29,   0x8000b,    0x200881,   0x288081,   0x820db,    0x52e23,    0x2888f,
            0xaaa8b,    0x28838d,   0x20000,    0xa93ab,    0x2c18d,    0x2ba89,    0x183600,   0xea00b,    0x1a4ea0,   0xf8585
        };
        set<Kmer<>::Type> okValues (values, values + ARRAY_SIZE(values));

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

