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

#include <gatb/kmer/impl/DSKAlgorithm.hpp>

#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/collections/impl/Product.hpp>
#include <gatb/tools/collections/impl/CollectionFile.hpp>

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

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestDSK : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestDSK);

        CPPUNIT_TEST_GATB (DSK_check1);
        CPPUNIT_TEST_GATB (DSK_check2);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    () {}
    void tearDown () {}

    /********************************************************************************/
    void DSK_check1_aux (const char* sequences[], size_t nbSequences, size_t kmerSize, size_t nks, size_t checkNbSolids)
    {
        /** We create a product instance. */
        Product<CollectionFile> product ("test");

        /** We create a DSK instance. */
        DSKAlgorithm<NativeInt64> dsk (product, new BankStrings (sequences, nbSequences), kmerSize, nks);

        /** We launch DSK. */
        dsk.execute();

        CPPUNIT_ASSERT (checkNbSolids == dsk.getInfo()->getInt("solid kmers nb"));
    }

    /********************************************************************************/
    void DSK_check1 ()
    {
        const char* s1 = "GATCCTCCCCAGGCCCCTACACCCAAT" ;

        const char* seqs1[] =  { s1 };
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 27, 1, 1);
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 26, 1, 2);
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 27, 2, 0);
        DSK_check1_aux (seqs1, ARRAY_SIZE(seqs1), 26, 2, 0);

        const char* seqs2[] =  { s1, s1 };
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 27, 1, 1);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 26, 1, 2);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 27, 2, 1);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 26, 2, 2);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 27, 3, 0);
        DSK_check1_aux (seqs2, ARRAY_SIZE(seqs2), 26, 3, 0);

        const char* seqs3[] =  { s1, s1, s1};
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 1, 1);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 1, 2);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 2, 1);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 2, 2);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 3, 1);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 3, 2);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 27, 4, 0);
        DSK_check1_aux (seqs3, ARRAY_SIZE(seqs3), 26, 4, 0);

        const char* seqs4[] = {

            "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
            "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
            "ACTTAGATGTAAGATTTCGAAGACTTGGATGTAAACAACAAATAAGATAATAACCATAAAAATAGAAATG"
            "AACGATATTAAAATTAAAAAATACGAAAAAACTAACACGTATTGTGTCCAATAAATTCGATTTGATAATT"
            "AGGTAACAATTTAACGTTAAAACCTATTCTTTTATTATCCGAAAATCCGTCGTGGAATTTGTATTAGCTT"
            "TTTTTCTACATTACCCGTTTGCGAGACAGGTGGGGTCAGACGTAGACGTAGTCTCTGGAGTCAAGACGAA"
            "ATTTTACATTTCACAATTTCCTATAGGCCGAGCAAAATTTATTAAGAACCCACAGGCATCATTACGTTTT"
            "CTTGCACAGAAGACTTCACGCTGAAGTCATTGGGCTATATTTCAACGAGACGTCTGTTGGTTTATAAAGG"
            "GCTATATTTATACAAGAATAGGAGTATGGCAGTATGCTAGGCTGGTATGTAGTACGTATACCTCCTAAGC"
            "CGAAAGGCAGTAAGTGACGATGTAATAGTTTTGAGGAAAATTACTTTTTCTGAATAATATTTTTATTTTT"
            "GTTTGCATTTTGTTAAAATTATTTACTAAATTAATGATTCTCATATGTTTTTTCATAGATTTGATGAACT"
            "ACTGTACCATCTGATTAGCGCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCA"
            "CACAACATACGAGCCGGAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTACTCACATAATTG"
            "CGTGCGCTCACTGCCCGCTTTCCAGTCGGAAACCTGTCGTGCCAGCTGCATTATTGATCCGCCACGCCCC"
            "GGGGAAAGT",

            "GTCTTCATTCAGCTGTTCTCATGATAACTAGTAATTCCTTGCTAACAATTTTTACTGAGTAGCAACCAAT"
            "TAATGTTGCCAGAATTTCATAATTGAATTTGAATTTTTTATTTTTTCCTTGATTATGCTTCAAACTCTAT"
            "GTAGTTATTTAGAGTCAATAATATTAAAGCAATCTTAATATTAACTCATTTATTTCTGATTGGCCATATT"
            "TATTTAATTCTCAACAATAATAATGATAAGTATAATAATATATTTAACTTAATAACATTTTAATCATTTT"
            "ATTTTTGTTTGTTGTGATTTTTGGACGTTGTGGTAAATAAGAAGTTTTAAGCTTATATTAATATGTTTTA"
            "CTTTTTATTTCTTAATACGAATTTAATTACCTACCCATTATATTAAGTATATGTTTTGGAATTCTTTCTG"
            "TAAAAATGTGTTTTAAATATTTTACACTTAATTATGTAGGTACCTATACATTTTTAAACTTATCGTATAA"
            "TTCTTTTAATGGTTAAATCATACAAATTAATGTGTAGAGAATAGTTTTTATAAGACTCGTTGTCAATACG"
            "TACGCATAATATAAAAAAACTGACATGTTTTAGTAAGTCGTTTTGATGCATAATAGGATTTTTACCTTTT"
            "AAAGTCTCAAGTTTTCATACAGTGGTACCTCTATATAGAACACGTTAGGCTTTACGGGGTCATTATTTCT"
            "GTTCCGATATTTTTAATGGCATAAAACTATAAACAATAACCGGTATGTATAAATGGTAC",

            "ACCTGAAGCAGTAGTTCATCATATTGCGACTGCAGAATCGATGATAAAGTGGCTTTTAGATCTAAAAGCC"
            "AATACAAAACTGAAGGAATTTGATTTGATGGATTTTAATTTTGAAAATGGATTATGATTGTCGATTGATT"
            "AACAAGTTTACTAGGTTTGAATAGAGGTGATTCTTAATATTTCAAATATTTGAATGTCATGATGAATATT"
            "ATAATTTATAATTAAAAAATATCATATTTTATTCATGGATATCAAAGCTGAAAAAATAGATATTCAAAAT"
            "CGCCTTTATAATAACCTATCATAAACTAATTAATCAATTAAATTCAGTTTTAAAAATTTAAATCCGACAA"
            "ATAAAATTCCTTCAGCTCTGTCTGGGATTTTGGTCGAAAAATTTTAAATCGAAAAAAGTTTATCTTATTC"
            "ATAATATCATTGCCAATGATATTAAAATTAATTAACAACGAATACAAATAACGTCCGACCTGTATATTGC"
            "GGGCCAACTGTTTTTATAGGAAATGTTGACCGAAAACTATTACAGATTAGATGTGTGTGTGTTTACCCTG"
            "TACAAAAATACAAGTACTATTACAACACATCATTATGTTAAATTGCCTCTATATTAATTTCTTTAAAACA"
            "CGACCAACTGCACATTAAAGTAGTTTATTTAGTACTACAGTAGATTAAATTCATTTTTGACGAAAAATTG"
            "CATTTGAAAATGGCCATTGTGTGTATAAATATTGTATACTAATATAACTCTAAATAAAGGTTTCCAGTAC"
            "CAAAGAACCAAATTTTTAATTACAACCAAAATAACTAAATCGTATTCTTTGTTAAATAGTTAAGTTTTTC"
            "GCCGATTGCTGTGCTTGACAGTCTCCTCAATTCAGAATTTCATGTAAAATAAAAATAGCGTACATATAAT"
            "GGATTGCTGTGGCATTTGGTTTGATTAATCCCAAATATTGATTCCAAATATCTATTAGCCTATTGTACCC"
            "CGGAGTACCG"
        } ;

        /** The following solid kmers numbers are computed with the original minia. */
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 1, 2540);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 2, 151);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 3, 18);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 4, 3);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 5, 2);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 9, 6, 0);

        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 11, 1, 2667);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 11, 2, 41);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 11, 3, 0);

        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 13, 1, 2690);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 13, 2, 12);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 13, 3, 0);

        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 15, 1, 2691);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 15, 2, 5);
        DSK_check1_aux (seqs4, ARRAY_SIZE(seqs4), 15, 3, 0);
    }

    /********************************************************************************/
    void DSK_check2 ()
    {
        size_t kmerSize = 31;
        size_t nks      = 1;

        const char* s1 = "GATCGATTCTTAGCACGTCCCCCCCTACACCCAAT" ;

        /** We create a product instance. */
        Product<CollectionFile> product ("test");

        /** We create a DSK instance. */
        DSKAlgorithm<NativeInt64> dsk (product, new BankStrings (s1, 0), kmerSize, nks);

        /** We launch DSK. */
        dsk.execute();

        /** We iterate the solid kmers. */
        Iterator<NativeInt64>* iter = dsk.getSolidKmers()->iterator();
        LOCAL (iter);

        /** The following values have been computed with the original DSK.
         *  Note: we use a set and check that the solid iterator items are in this trustful set.
         *  -> we have to do this because we are not sure about the order of the iterated items.
         */
        set<NativeInt64> okValues;
        okValues.insert (0x1CA68D1E55561150);
        okValues.insert (0x09CA68D1E5556115);
        okValues.insert (0x2729A34795558454);
        okValues.insert (0x32729A3479555845);
        okValues.insert (0x0AFEE3FFF1ED8309);

        set<NativeInt64> checkValues;

        size_t idx=0;
        for (iter->first(); !iter->isDone(); iter->next(), idx++)
        {
            set<NativeInt64>::iterator lookup = okValues.find (iter->item());
            CPPUNIT_ASSERT (lookup != okValues.end());

            checkValues.insert (iter->item());
        }

        CPPUNIT_ASSERT (checkValues.size() == okValues.size());
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestDSK);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestDSK);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

