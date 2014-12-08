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
#include <gatb/tools/misc/api/Data.hpp>
#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/bank/api/Sequence.hpp>
#include <gatb/bank/impl/Alphabet.hpp>
#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankRandom.hpp>
#include <gatb/bank/impl/BankStrings.hpp>
#include <gatb/tools/collections/impl/IterableHelpers.hpp>

#include <gatb/system/api/Exception.hpp>

#include <iostream>

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

extern std::string DBPATH (const string& a);

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
        CPPUNIT_TEST_GATB (kmer_minimizer);
        CPPUNIT_TEST_GATB (kmer_minimizer2);
        CPPUNIT_TEST_GATB (kmer_badchar);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    template<size_t span> class Check
    {
    public:

        typedef typename Kmer<span>::ModelCanonical Model;
        typedef typename Kmer<span>::Type           Type;

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
        template <class ModelType>
        static void kmer_checkCompute_aux (const char* seq, long* check, size_t len)
        {
            typename ModelType::Kmer kmer;

            /** We declare a kmer model with a given span size. */
            ModelType model (3);

            /** We compute the kmer for a given sequence */
            kmer = model.codeSeed (seq, Data::ASCII);
            CPPUNIT_ASSERT (kmer.value() == check[0]);

            /** We compute some of the next kmers. */
            kmer = model.codeSeedRight (kmer, seq[3], Data::ASCII);

            CPPUNIT_ASSERT (kmer.value() == check[1]);

            kmer = model.codeSeedRight (kmer, seq[4], Data::ASCII);
            CPPUNIT_ASSERT (kmer.value() == check[2]);

            kmer = model.codeSeedRight (kmer, seq[5], Data::ASCII);
            CPPUNIT_ASSERT (kmer.value() == check[3]);
        }

        /** */
        static void kmer_checkCompute ()
        {
            const char* seq = "CATTGATAGTGG";
            long directKmers [] = {18, 10, 43, 44, 50,  8, 35, 14, 59, 47};
            long reverseKmers[] = {11,  2, 16, 36,  9, 34, 24,  6, 17, 20};

            kmer_checkCompute_aux <typename Kmer<span>::ModelDirect> (seq, directKmers,  ARRAY_SIZE(directKmers));
        }

        /********************************************************************************/
        template<class ModelType>
        static void kmer_checkIterator_aux (const char* seq, long* kmersTable, size_t lenTable)
        {
            /** We declare a kmer model with a given span size. */
            ModelType model (3);

            /** We declare an iterator. */
            typename ModelType::Iterator it (model);

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
                CPPUNIT_ASSERT (it->value() == kmersTable[idx++]);
            }
            /** We check we found the correct number of kmers. */
            CPPUNIT_ASSERT (idx == lenTable);
        }

        /********************************************************************************/
        static void kmer_checkIterator ()
        {
            const char* seq = "CATTGATAGTGG";

             long checkDirect []  = {18, 10, 43, 44, 50,  8, 35, 14, 59, 47};
             kmer_checkIterator_aux <typename Kmer<span>::ModelDirect> (seq, checkDirect, ARRAY_SIZE(checkDirect));

             // long checkReverse [] = {11,  2, 16, 36,  9, 34, 24,  6, 17, 20};
             // kmer_checkIterator_aux (model, seq, KMER_REVCOMP, checkReverse, ARRAY_SIZE(checkReverse));

            long checkBoth []    = {11,  2, 16, 36,  9,  8, 24,  6, 17, 20};
            kmer_checkIterator_aux <typename Kmer<span>::ModelCanonical> (seq, checkBoth, ARRAY_SIZE(checkBoth));
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

        Kmer<>::ModelCanonical model (5);

        Kmer<>::Type check[] = {0x61, 0x187, 0x21c, 0x72, 0x1c9, 0x1c9, 0x9c, 0x9c, 0x127, 0x49, 0xb8};

        vector<Kmer<>::ModelCanonical::Kmer> kmers;

        model.build (s1.getData(), kmers);

        CPPUNIT_ASSERT (kmers.size() == 11);

        size_t i=0;
        for (i=0; i<kmers.size(); i++)  {  CPPUNIT_ASSERT (kmers[i].value() == check[i]);  }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(check));

        for (i=0; i<ARRAY_SIZE(check); i++)  {  CPPUNIT_ASSERT (model.getKmer(s1.getData(), i).value() == check[i]);  }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(check));

        for (i=0; i<ARRAY_SIZE(check); i++)  {  CPPUNIT_ASSERT (model.getKmer (Data(buf), i).value() == check[i]);  }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(check));

        Kmer<>::ModelCanonical::Kmer kmer = model.getKmer (Data(buf));
        CPPUNIT_ASSERT (kmer.value() == check[0]);
    }

    /********************************************************************************/
    template<size_t span>
    struct kmer_minimizer_fct
    {
        typedef typename Kmer<span>::ModelDirect ModelDirect;
        typedef typename Kmer<span>::template ModelMinimizer<ModelDirect> ModelMinimizer;

        ModelMinimizer& model;
        size_t& nbKmers;

        kmer_minimizer_fct (ModelMinimizer& model, size_t& nbKmers)  : model(model), nbKmers(nbKmers) {}

        void operator() (const typename ModelMinimizer::Kmer& kmers, size_t idx)
        {
            string currentKmer = model.toString(kmers.value());

            typename Kmer<span>::Type currentMini (~0);
            for (size_t i=0; i<model.getKmerSize() - model.getMmersModel().getKmerSize() + 1; i++)
            {
                typename ModelDirect::Kmer tmp = model.getMmersModel().codeSeed (
                    currentKmer.substr(i,model.getMmersModel().getKmerSize()).data(),
                    Data::ASCII
                );

                if (tmp.value() < currentMini)  { currentMini = tmp.value(); }
            }

            CPPUNIT_ASSERT (currentMini == kmers.minimizer().value() );

            nbKmers++;
        }
    };

    /** */
    template<size_t span, class ModelType>
    void kmer_minimizer_aux2 (IBank& bank, size_t kmerSize, size_t miniSize)
    {
        typename Kmer<span>::template ModelMinimizer <ModelType> minimizerModel (kmerSize, miniSize);

        size_t nbKmers = 0;
        // We loop over sequences.
        Iterator<Sequence>* itSeq = bank.iterator();  LOCAL (itSeq);

        for (itSeq->first(); !itSeq->isDone(); itSeq->next())
        {
            minimizerModel.iterate ((*itSeq)->getData(), kmer_minimizer_fct<span> (minimizerModel, nbKmers));
        }
        CPPUNIT_ASSERT (nbKmers > 0);
    }

    /** */
    void kmer_minimizer ()
    {
        /** NOTE: sequences should be big enough to get kmers up to 128. */
        vector<IBank*> banks;
        banks.push_back (new BankStrings ("ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAGAT"
                                          "TCCCCAGGCCCCTACACCCAATGTGGAACCGGGGTCCCGAATGAAAATGCTGCTGTTCCCTGGAGGTGTTCT", NULL));
        banks.push_back (Bank::open (DBPATH("reads1.fa")));
        banks.push_back (new BankRandom (500, 200));

        size_t   kmerSizes[] = { 12,31, 33,63, 65,95, 97,127 };
        size_t   miniSizes[] = {  5,  9, 11 };

        try
        {
            for (size_t b=0; b<banks.size(); b++)
            {
                IBank* bank = banks[b];   LOCAL(bank);

                for (size_t i=0; i<ARRAY_SIZE(kmerSizes); i++)
                {
                    size_t kmerSize = kmerSizes[i];

                    for (size_t j=0; j<ARRAY_SIZE(miniSizes); j++)
                    {
                             if (kmerSize < KSIZE_1)  {  kmer_minimizer_aux2<KSIZE_1, Kmer<KSIZE_1>::ModelDirect> (*bank, kmerSize, miniSizes[j]);  }
                        else if (kmerSize < KSIZE_2)  {  kmer_minimizer_aux2<KSIZE_2, Kmer<KSIZE_2>::ModelDirect> (*bank, kmerSize, miniSizes[j]);  }
                        else if (kmerSize < KSIZE_3)  {  kmer_minimizer_aux2<KSIZE_3, Kmer<KSIZE_3>::ModelDirect> (*bank, kmerSize, miniSizes[j]);  }
                    }
                }
            }
        }
        catch (gatb::core::system::Exception& e)
        {
            cout << e.getMessage() << endl;
        }
    }

    /** */
    struct kmer_minimizer2_info
    {
        const char* kmer;
        const char* minimizer;
        size_t      position;
        bool        changed;
    };

    void kmer_minimizer2 ()
    {
        const char* seq = "ATGTCTGAAGTGACCTAACATTGCA";
        size_t sequenceSize = strlen(seq);

        size_t kmerSize = 15;
        size_t mmerSize = 7;

        typedef Kmer<KSIZE_1>::ModelDirect                  ModelDirect;
        typedef Kmer<KSIZE_1>::ModelMinimizer<ModelDirect>  ModelMinimizer;

        ModelMinimizer model (kmerSize, mmerSize);

        const ModelDirect& modelMini = model.getMmersModel();

        kmer_minimizer2_info table[] =
        {
            {"ATGTCTGAAGTGACC", "AAGTGAC", 7, true },
            {"TGTCTGAAGTGACCT", "AAGTGAC", 6, false},
            {"GTCTGAAGTGACCTA", "AAGTGAC", 5, false},
            {"TCTGAAGTGACCTAA", "AAGTGAC", 4, false},
            {"CTGAAGTGACCTAAC", "AAGTGAC", 3, false},
            {"TGAAGTGACCTAACA", "AAGTGAC", 2, false},
            {"GAAGTGACCTAACAT", "AAGTGAC", 1, false},
            {"AAGTGACCTAACATT", "AAGTGAC", 0, false},
            {"AGTGACCTAACATTG", "AACATTG", 8, true },
            {"GTGACCTAACATTGC", "AACATTG", 7, false},
            {"TGACCTAACATTGCA", "AACATTG", 6, false}
          };

        size_t idx = 0;

        #define CHECK(idx) \
            CPPUNIT_ASSERT (  \
                   model.toString (kmer.value()) == table[idx].kmer \
                && modelMini.toString(kmer.minimizer().value()) == table[idx].minimizer \
                && kmer.position()   == table[idx].position  \
                && kmer.hasChanged() == table[idx].changed)

        ModelMinimizer::Kmer kmer = model.codeSeed (seq, Data::ASCII);
        CHECK (idx);  idx++;
        seq += kmerSize;

        for ( ; idx < ARRAY_SIZE(table); idx++)
        {
            kmer = model.codeSeedRight (kmer, *(seq++), Data::ASCII);
            CHECK (idx);
        }
    }

    /********************************************************************************/

    typedef Kmer<KSIZE_1>::ModelDirect  ModelDirect;

    struct kmer_badchar_info
    {
        const char* kmer;
        bool        valid;
    };

    struct kmer_badchar_functor
    {
        const ModelDirect& model;
        kmer_badchar_info* table;
        size_t             length;

        kmer_badchar_functor (const ModelDirect& model, kmer_badchar_info* table, size_t length) : model(model), table(table), length(length) {}

        void operator() (const ModelDirect::Kmer& kmer, size_t idx)
        {
            CPPUNIT_ASSERT (kmer.isValid() == table[idx].valid);

            if (kmer.isValid())
            {
                CPPUNIT_ASSERT (model.toString(kmer.value()) == table[idx].kmer);
            }
            else
            {
                // Bad nucleotides should have replaced N by G
                string modif (table[idx].kmer);
                for (size_t i=0; i<modif.size(); i++)  { if (modif[i]=='N') { modif[i]='G'; }}
                CPPUNIT_ASSERT (model.toString(kmer.value()) == modif);
            }
        }
    };

    /** */
    void kmer_badchar (void)
    {
        typedef Kmer<KSIZE_1>::ModelDirect  ModelDirect;

        size_t kmerSize = 11;
        ModelDirect model (kmerSize);

        const char* seq = "ACGNCNTGCTAGCTATTTAGCTTTAGANAGTAGATGACGCNC";

        kmer_badchar_info info[] =
        {
            {"ACGNCNTGCTA", false}, {"CGNCNTGCTAG", false}, {"GNCNTGCTAGC", false}, {"NCNTGCTAGCT", false},
            {"CNTGCTAGCTA", false}, {"NTGCTAGCTAT", false}, {"TGCTAGCTATT", true }, {"GCTAGCTATTT", true },
            {"CTAGCTATTTA", true }, {"TAGCTATTTAG", true }, {"AGCTATTTAGC", true }, {"GCTATTTAGCT", true },
            {"CTATTTAGCTT", true }, {"TATTTAGCTTT", true }, {"ATTTAGCTTTA", true }, {"TTTAGCTTTAG", true },
            {"TTAGCTTTAGA", true }, {"TAGCTTTAGAN", false}, {"AGCTTTAGANA", false}, {"GCTTTAGANAG", false},
            {"CTTTAGANAGT", false}, {"TTTAGANAGTA", false}, {"TTAGANAGTAG", false}, {"TAGANAGTAGA", false},
            {"AGANAGTAGAT", false}, {"GANAGTAGATG", false}, {"ANAGTAGATGA", false}, {"NAGTAGATGAC", false},
            {"AGTAGATGACG", true }, {"GTAGATGACGC", true }, {"TAGATGACGCN", false}, {"AGATGACGCNC", false}
        };

        Data data (Data::ASCII);
        data.set ((char*)seq, strlen(seq));

        kmer_badchar_functor fct (model, info, ARRAY_SIZE(info));

        model.iterate (data, fct);
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestKmer);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestKmer);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

