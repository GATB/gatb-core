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
#include <gatb/tools/misc/api/Data.hpp>
#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/bank/api/Sequence.hpp>
#include <gatb/bank/impl/Alphabet.hpp>
#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>

#include <iostream>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestKmer : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestKmer);

        CPPUNIT_TEST_GATB (kmer_checkInfo);
        CPPUNIT_TEST_GATB (kmer_checkCompute);
        CPPUNIT_TEST_GATB (kmer_checkIterator);
        CPPUNIT_TEST_GATB (kmer_build);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    template<size_t span> class Check
    {
    public:

        typedef typename Kmer<span>::Model Model;
        typedef typename Kmer<span>::Type  Type;

        /********************************************************************************/
        static void kmer_checkInfo ()
        {
            size_t aSpan = 27;

            /** We declare a kmer model with a given span size. */
            Model model (aSpan);

            /** We check some information. */
            CPPUNIT_ASSERT (model.getKmerSize()    == aSpan);
            CPPUNIT_ASSERT (model.getMemorySize()  == sizeof(Type));
        }

        /********************************************************************************/
        static void kmer_checkCompute_aux (const char* seq, long* check, size_t len, KmerMode mode)
        {
            Type kmer;

            /** We declare a kmer model with a given span size. */
            Model model (3, mode);

            /** We compute the kmer for a given sequence */
            kmer = model.codeSeed (seq, Data::ASCII);
            CPPUNIT_ASSERT (kmer == check[0]);

            /** We compute some of the next kmers. */
            kmer = model.codeSeedRight (kmer, seq[3], Data::ASCII);
            CPPUNIT_ASSERT (kmer == check[1]);

            kmer = model.codeSeedRight (kmer, seq[4], Data::ASCII);
            CPPUNIT_ASSERT (kmer == check[2]);

            kmer = model.codeSeedRight (kmer, seq[5], Data::ASCII);
            CPPUNIT_ASSERT (kmer == check[3]);
        }

        /** */
        static void kmer_checkCompute ()
        {
            const char* seq = "CATTGATAGTGG";
            long directKmers [] = {18, 10, 43, 44, 50,  8, 35, 14, 59, 47};
            long reverseKmers[] = {11,  2, 16, 36,  9, 34, 24,  6, 17, 20};

            kmer_checkCompute_aux (seq, directKmers,  ARRAY_SIZE(directKmers),  KMER_DIRECT);
        }

        /********************************************************************************/
        static void kmer_checkIterator_aux (Model& model, const char* seq, KmerMode mode, long* kmersTable, size_t lenTable)
        {
            /** We declare an iterator. */
            typename Model::Iterator it (model);

            /** We set the data from which we want to extract kmers. */
            Data data (Data::ASCII);
            data.set ((char*)seq, strlen(seq));

            /** We feed the iterator with the sequence data content. */
            it.setData (data);
            /** We iterate the kmers. */
            size_t idx=0;
            for (it.first(); !it.isDone(); it.next())
            {
                /** We check that the iterated kmer is the good one. */
                CPPUNIT_ASSERT (*it == kmersTable[idx++]);
            }
            /** We check we found the correct number of kmers. */
            CPPUNIT_ASSERT (idx == lenTable);
        }

        /********************************************************************************/
        static void kmer_checkIterator ()
        {
            /** We declare a kmer model with a given span size. */
            typename Kmer<span>::Model model (3);

            const char* seq = "CATTGATAGTGG";

            // long checkDirect []  = {18, 10, 43, 44, 50,  8, 35, 14, 59, 47};
            // kmer_checkIterator_aux (model, seq, KMER_DIRECT, checkDirect, ARRAY_SIZE(checkDirect));
            //
            // long checkReverse [] = {11,  2, 16, 36,  9, 34, 24,  6, 17, 20};
            // kmer_checkIterator_aux (model, seq, KMER_REVCOMP, checkReverse, ARRAY_SIZE(checkReverse));

            long checkBoth []    = {11,  2, 16, 36,  9,  8, 24,  6, 17, 20};
            kmer_checkIterator_aux (model, seq, KMER_MINIMUM, checkBoth, ARRAY_SIZE(checkBoth));
        }
    };

    /********************************************************************************/
    void kmer_checkInfo ()
    {
        /** We check with the native 64 bits type. */
        Check<32>::kmer_checkInfo();

        /** We check with the native 128 bits type. */
        Check<64>::kmer_checkInfo();

        /** We check with the LargeInt type. */
        Check<96>::kmer_checkInfo();
    }

    /********************************************************************************/
    void kmer_checkCompute ()
    {
        /** We check with the native 64 bits type. */
        Check<32>::kmer_checkCompute();

        /** We check with the native 128 bits type. */
        Check<64>::kmer_checkCompute();

        /** We check with the LargeInt type. */
        Check<96>::kmer_checkCompute();
    }

    /********************************************************************************/
    void kmer_checkIterator ()
    {
        /** We check with the native 64 bits type. */
        Check<32>::kmer_checkIterator();

        /** We check with the native 128 bits type. */
        Check<64>::kmer_checkIterator();

        /** We check with the LargeInt type. */
        Check<96>::kmer_checkIterator();
    }

    /********************************************************************************/
    void kmer_build ()
    {
        char* buf = (char*)"ACTACGATCGATGTA";

        Sequence s1 (buf);

        Kmer<>::Model model (5);

        Kmer<>::Type check[] = {0x61, 0x187, 0x21c, 0x72, 0x1c9, 0x1c9, 0x9c, 0x9c, 0x127, 0x49, 0xb8};

        vector<Kmer<>::Type> kmers;

        model.build (s1.getData(), kmers);
        CPPUNIT_ASSERT (kmers.size() == 11);

        size_t i=0;
        for (i=0; i<kmers.size(); i++)  {  CPPUNIT_ASSERT (kmers[i] == check[i]);  }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(check));

        for (i=0; i<ARRAY_SIZE(check); i++)  {  CPPUNIT_ASSERT (model.getKmer(s1.getData(), i).first == check[i]);  }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(check));

        for (i=0; i<ARRAY_SIZE(check); i++)  {  CPPUNIT_ASSERT (model.getKmer (Data(buf), i).first == check[i]);  }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(check));

        Kmer<>::Type kmer = model.getKmer (Data(buf)).first;
        CPPUNIT_ASSERT (kmer == check[0]);
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestKmer);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestKmer);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

