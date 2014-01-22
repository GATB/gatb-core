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

#include <gatb/bank/impl/Banks.hpp>

#include <gatb/bank/impl/BankRegistery.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/misc/api/Macros.hpp>

#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

#define ABS(a)  ((a)<0 ? -(a) : (a))

extern std::string DBPATH (const string& a);

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestBank : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestBank);

        CPPUNIT_TEST_GATB (bank_checkSample1);
        CPPUNIT_TEST_GATB (bank_checkSample2);
        CPPUNIT_TEST_GATB (bank_checkComments);
        CPPUNIT_TEST_GATB (bank_checkBadUri);
        CPPUNIT_TEST_GATB (bank_checkSize);
        CPPUNIT_TEST_GATB (bank_checkNbFiles);
        CPPUNIT_TEST_GATB (bank_checkEstimateNbSequences);
        CPPUNIT_TEST_GATB (bank_checkProgress);
        CPPUNIT_TEST_GATB (bank_checkConvertBinary);
        CPPUNIT_TEST_GATB (bank_checkRegistery1);
        CPPUNIT_TEST_GATB (bank_checkRegistery2);
        CPPUNIT_TEST_GATB (bank_strings1);
        CPPUNIT_TEST_GATB (bank_sequence);
        CPPUNIT_TEST_GATB (bank_splitter_1);
        CPPUNIT_TEST_GATB (bank_random_1);
        CPPUNIT_TEST_GATB (bank_composite);
        CPPUNIT_TEST_GATB (bank_album1);
        CPPUNIT_TEST_GATB (bank_album2);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/
    void bank_checkSample1_aux (const string& filename, BankFasta::Iterator::CommentMode_e mode)
    {
        const char* text = "ARNDCQEGHILKMFPSTWYV";

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b (filename);

        /** We create an iterator over this bank. */
        BankFasta::Iterator it (b, mode);

        size_t i=0;

        /** We loop over sequences. */
        for (it.first(); !it.isDone(); it.next(), i++)
        {
            /** We prepare the sequence string to be matched according to the provided mode. */
            char buffer[32];
            snprintf (buffer, sizeof(buffer), "seq%ld%s", i+1,  (mode == BankFasta::Iterator::FULL ? " generic" : ""));

            /** We check that we got a comment and that it is like 'seqX' where X is the
             * index from 1 to 20. */
            CPPUNIT_ASSERT (it->getComment().empty() == false);
            CPPUNIT_ASSERT (it->getComment().compare(buffer) == 0);

            /** We check that the data size is 20. */
            CPPUNIT_ASSERT (it->getDataSize() == 20);

            /** We check that the retrieved sequence is compliant with the sample database. */
            for (size_t k=0; k<20; k++)  {  CPPUNIT_ASSERT (it->getData()[k] == text[ (i+k) % 20]);  }
        }

        /** We check that we read exactly 20 sequences. */
        CPPUNIT_ASSERT (i == 20);
    }

    /********************************************************************************/
    /** \brief check Bank class information on a specific bank.
     *
     * In this test, we use a specific 'sample' FASTA bank:
     *      - it holds 20 sequences
     *      - sequences comment are like seq1, seq2, ... seq20
     *      - data of seq1 is ARNDCQEGHILKMFPSTWYV
     *      - data of seq2 is RNDCQEGHILKMFPSTWYVA
     *      - other sequences continue a cyclic left rotation of one letter
     *
     * We iterate this sample bank and check that iterated gatb::core::bank::Sequence
     * instances are coherent with the specification of the sample bank.
     *
     * The test is done both with the uncompressed and compressed banks.
     *
     * Test of \ref gatb::core::bank::impl::Bank                \n
     * Test of \ref gatb::core::bank::impl::Bank::Iterator      \n
     * Test of \ref gatb::core::bank::Sequence                  \n
     * Test of \ref gatb::core::bank::Sequence::getComment()    \n
     * Test of \ref gatb::core::bank::Sequence::getDataSize()   \n
     * Test of \ref gatb::core::bank::Sequence::getData()       \n
     */
    void bank_checkSample1 ()
    {
        /** We launch the test with uncompressed sample bank. */
        bank_checkSample1_aux (DBPATH("sample1.fa"), BankFasta::Iterator::IDONLY);
        bank_checkSample1_aux (DBPATH("sample1.fa"), BankFasta::Iterator::FULL);

        /** We launch the test with compressed sample bank. */
        bank_checkSample1_aux (DBPATH("sample1.fa.gz"), BankFasta::Iterator::IDONLY);
        bank_checkSample1_aux (DBPATH("sample1.fa.gz"), BankFasta::Iterator::FULL);
    }

    /********************************************************************************/
    void bank_checkSample2_aux (const string& filename)
    {
        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b (filename);

        /** We create an iterator over this bank with only id comments. */
        BankFasta::Iterator it (b, BankFasta::Iterator::IDONLY);

        size_t nbSeq=0;

        /** We loop over sequences. */
        for (it.first(); !it.isDone(); it.next(), nbSeq++)
        {
            char buffer[32];
            snprintf (buffer, sizeof(buffer), "seq%ld", nbSeq+1);

            /** We check that we got a comment and that it is like 'seqX' where X is the
             * index from 1 to 20. */
            CPPUNIT_ASSERT (it->getComment().empty() == false);
            CPPUNIT_ASSERT (it->getComment().compare (buffer) == 0);

            /** We check that the data size is 0. */
            CPPUNIT_ASSERT (it->getDataSize() == 0);
        }

        /** We check that we read exactly 20 sequences. */
        CPPUNIT_ASSERT (nbSeq == 20);
    }

    /********************************************************************************/
    /** \brief check Bank class information on a bank with comments but without data.
     *
     * Test of \ref gatb::core::bank::impl::Bank                \n
     * Test of \ref gatb::core::bank::impl::Bank::Iterator      \n
     * Test of \ref gatb::core::bank::Sequence                  \n
     * Test of \ref gatb::core::bank::Sequence::getComment()    \n
     * Test of \ref gatb::core::bank::Sequence::getDataSize()   \n
     */
    void bank_checkSample2 ()
    {
        /** We launch the test with uncompressed sample bank. */
        bank_checkSample2_aux (DBPATH("sample2.fa"));
    }

    /********************************************************************************/
    void bank_checkComments_aux (const string& filename)
    {
        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b1 (filename);

        /** We iterate without comments. */
        for (BankFasta::Iterator it (b1, BankFasta::Iterator::NONE); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty());
        }

        /** We iterate with only id comments. */
        for (BankFasta::Iterator it (b1, BankFasta::Iterator::IDONLY); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty() == false);
            CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") == 0);
        }

        /** We iterate with only full comments. */
        for (BankFasta::Iterator it (b1, BankFasta::Iterator::FULL); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty() == false);
            CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") != 0);
        }

        /** We iterate with default value (should be FULL). */
        for (BankFasta::Iterator it (b1); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty()==false);
            CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") != 0);
        }
    }

    /********************************************************************************/
    /** \brief check how Bank::Iterator class manages sequence comments
     *
     * A database is iterated 3 times:
     *      1) without retrieving sequence comments
     *      2) with retrieving only sequence id
     *      3) with retrieving full sequence comments
     *
     * We check in each case the coherence of the retrieved sequence comments.
     *
     * Test of \ref gatb::core::bank::impl::Bank                \n
     * Test of \ref gatb::core::bank::impl::Bank::Iterator      \n
     * Test of \ref gatb::core::bank::Sequence                  \n
     * Test of \ref gatb::core::bank::Sequence::getComment()    \n
     */
    void bank_checkComments ()
    {
        /** We launch the test with uncompressed bank. */
        bank_checkComments_aux (DBPATH("query.fa"));

        /** We launch the test with compressed bank. */
        bank_checkComments_aux (DBPATH("query.fa.gz"));
    }

    /********************************************************************************/
    /** \brief ok/ko bank uri check
     *
     * This test creates Bank instance with bad and good filenames:
     *      - bad test:   check we get a gatb::core::system::Exception exception
     *      - good test:  check we get no exception
     *
     * Test of \ref gatb::core::bank::impl::Bank
     */
    void bank_checkBadUri ()
    {
        BankFasta bankKO ("_dummy_bank_name_");
        BankFasta bankOK (DBPATH("sample1.fa"));

        /** We check that we can't get a valid iterator on it. */
        CPPUNIT_ASSERT_THROW (BankFasta::Iterator it(bankKO), gatb::core::system::Exception);

        /** We check that the bank does exist. */
        CPPUNIT_ASSERT_NO_THROW (BankFasta::Iterator it(bankOK));
    }

    /********************************************************************************/
    /** \brief Check the size of a bank file
     *
     * We just check that we can retrieve the correct (and known) size of one sample bank.
     *
     * Test of \ref gatb::core::bank::impl::Bank                \n
     * Test of \ref gatb::core::bank::impl::Bank::getSize()     \n
     */
    void bank_checkSize ()
    {
        /** We declare a Bank instance. */
        BankFasta b1 (DBPATH("sample1.fa"));

        /** We check the size of the bank. */
        CPPUNIT_ASSERT (b1.getSize() == 710);
    }

    /********************************************************************************/
    /** \brief Check the number of file names for creating a bank
     *
     * We check that we can't create a Bank object with 0 or too many filenames.
     * We check also that everything is ok between these two limits.
     *
     * Test of \ref gatb::core::bank::impl::Bank                    \n
     * Test of \ref gatb::core::bank::impl::Bank::getMaxNbFiles()   \n
     */
    void bank_checkNbFiles ()
    {
        vector<string> filenames;

        /** We check that we can't use a Bank without at least one filename. */
        CPPUNIT_ASSERT_THROW (BankFasta b (filenames), gatb::core::system::Exception);

        while (filenames.size() < BankFasta::getMaxNbFiles())
        {
            filenames.push_back(DBPATH("sample1.fa"));

            CPPUNIT_ASSERT_NO_THROW (BankFasta b (filenames));
        }

        /** We check that we can't use a Bank with with too many filenames. */
        filenames.push_back(DBPATH("sample1.fa"));
        CPPUNIT_ASSERT_THROW (BankFasta b (filenames), gatb::core::system::Exception);
    }

    /********************************************************************************/
    /** \brief Check that we can estimate the number of sequences in a bank
     *
     * Since we know the number of sequences in one sample bank, we can correctly
     * retrieve this number through the IBank API.
     *
     * The test is done with one simple sample, then with banks defined by several files.
     *
     * Test of \ref gatb::core::bank::impl::Bank                        \n
     * Test of \ref gatb::core::bank::impl::Bank::estimateNbSequences() \n
     * Test of \ref gatb::core::bank::impl::Bank::getMaxNbFiles()       \n
     */
    void bank_checkEstimateNbSequences ()
    {
        /** We declare a Bank instance. */
        BankFasta b1 (DBPATH("sample1.fa"));

        /** We check the estimation of sequences number. */
        u_int64_t estim1 = b1.estimateNbSequences();
        CPPUNIT_ASSERT (estim1 == 20);

        /** We build another bank holding several time the same bank. */
        for (size_t i=1; i<=BankFasta::getMaxNbFiles(); i++)
        {
            vector<string> filenames;
            for (size_t j=1; j<=i; j++)  { filenames.push_back(DBPATH("sample1.fa")); }

            BankFasta b2 (filenames);

            CPPUNIT_ASSERT (b2.estimateNbSequences() == i * estim1);
        }
    }

    /********************************************************************************/

    /** */
    struct ProgressFunctor : public IteratorListener
    {
        size_t nbCalls;
        size_t modulo;

        ProgressFunctor (size_t mod) : nbCalls(0), modulo(mod) {}

        void inc (u_int64_t current)   {  nbCalls++;  CPPUNIT_ASSERT (current % modulo == 0);  }
    };

    /** */
    void bank_checkProgress_aux (BankFasta& b, size_t modulo)
    {
        /** We create an iterator for tha provided bank. */
        Iterator<Sequence>* itSeq = b.iterator();
        LOCAL (itSeq);

        /** We create a progress iterator from a Sequence iterator provided by the bank. */
        SubjectIterator<Sequence> it (itSeq, modulo);

        /** We create some functor to be notified every N iteration and attach it to the iterator. */
        ProgressFunctor* fct = new ProgressFunctor (modulo);
        LOCAL (fct);
        it.addObserver (fct);

        /** We loop over the sequences; the functor should be called every N iterations. */
        for (it.first(); !it.isDone(); it.next())
        {
            /** NOTE ! we don't spoil the inner iteration loop with progression stuff. */
        }

        /** We check that our observer has been called as much as wanted.
         *  Note that we can rely here on the estimateNbSequences method since the bank is
         *  a known sample
         */
        CPPUNIT_ASSERT (modulo > 0);
        CPPUNIT_ASSERT (fct->nbCalls > 0);

        /** We have to check if we got the correct number of notifications. */
        CPPUNIT_ASSERT (fct->nbCalls == (b.estimateNbSequences() + modulo - 1) / modulo);

        /** We keep the number of functor calls. */
        size_t keepNbCalls = fct->nbCalls;

        /** We unsubscribe the functor from the iterator. */
        it.removeObserver (fct);

        /** We loop again over the sequences; the functor should not be called anymore. */
        for (it.first(); !it.isDone(); it.next())
        {
        }

        /** We check that the functor has not been called. */
        CPPUNIT_ASSERT (fct->nbCalls == keepNbCalls);
    }

    /** \brief Test the bank iterator with progression facilities
     *
     * This test creates "progress" iterator from a bank iterator and associate a listener to it.
     * We check that the listener is correctly called (knowing the properties of the used bank).
     *
     * Test of \ref gatb::core::bank::impl::Bank                                    \n
     * Test of \ref gatb::core::bank::impl::Bank::estimateNbSequences()             \n
     * Test of \ref gatb::core::tools::dp::impl::SubjectIterator                    \n
     * Test of \ref gatb::core::tools::dp::impl::SubjectIterator::addObserver       \n
     * Test of \ref gatb::core::tools::dp::impl::SubjectIterator::removeObserver    \n
     */
    void bank_checkProgress ()
    {
        string filename = DBPATH("sample1.fa");

        vector<string> filenames;

        size_t tableMod[] = { 1, 2, 3, 5, 8, 13};

        /** We build a bank holding x1 the same bank.*/
        while (filenames.size() < 1)   {  filenames.push_back (filename);  }
        BankFasta b1 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b1, tableMod[i]);  }

        /** We build a bank holding x2 the same bank.*/
        while (filenames.size() < 2)   {  filenames.push_back (filename);  }
        BankFasta b2 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b2, tableMod[i]);  }

        /** We build a bank holding x10 the same bank.*/
        while (filenames.size() < 10)   {  filenames.push_back (filename);  }
        BankFasta b10 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b10, tableMod[i]);  }

        /** We build a bank holding x20 the same bank.*/
        while (filenames.size() < 20)   {  filenames.push_back (filename);  }
        BankFasta b20 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b20, tableMod[i]);  }
    }

    /********************************************************************************/
    void bank_checkMultipleFiles_aux (const string& filename)
    {
        vector<string> filenames;

        /** A utility structure for this test. */
        struct Info  {  u_int32_t nbseq;  u_int64_t datasize;    Info() : nbseq(0), datasize(0) {}  };

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b1 (filename);

        /** We gather some information about this (single) bank. */
        Info i1;
        for (BankFasta::Iterator it (b1); !it.isDone(); it.next())   {  i1.nbseq++;   i1.datasize += it->getDataSize(); }

        /** We build a bank holding x2 the same bank.*/
        while (filenames.size() < 2)   {  filenames.push_back (filename);  }
        BankFasta b2 (filenames);

        /** We gather some information about this x2 bank. */
        Info i2;
        for (BankFasta::Iterator it (b2); !it.isDone(); it.next())   {  i2.nbseq++;   i2.datasize += it->getDataSize(); }

        /** We build a bank holding x4 the same bank.*/
        while (filenames.size() < 4)   {  filenames.push_back (filename);  }
        BankFasta b4 (filenames);

        /** We gather some information about this x4 bank. */
        Info i4;
        for (BankFasta::Iterator it (b4); !it.isDone(); it.next())   {  i4.nbseq++;   i4.datasize += it->getDataSize(); }

        CPPUNIT_ASSERT (i1.nbseq   * 2 == i2.nbseq);
        CPPUNIT_ASSERT (i1.datasize* 2 == i2.datasize);

        CPPUNIT_ASSERT (i1.nbseq   * 4 == i4.nbseq);
        CPPUNIT_ASSERT (i1.datasize* 4 == i4.datasize);

        CPPUNIT_ASSERT (i2.nbseq   * 2 == i4.nbseq);
        CPPUNIT_ASSERT (i2.datasize* 2 == i4.datasize);
    }

    /********************************************************************************/
    /** \brief Check that multiple files reading works well.
     *
     * This test takes a list of banks and creates an iterator on it. The idea is to
     * check that iterating N times the same bank with one iterator will provide N times
     * the information provided by iterating the bank only once.
     *
     * Test of \ref gatb::core::bank::impl::Bank                \n
     * Test of \ref gatb::core::bank::impl::Bank::Iterator      \n
     * Test of \ref gatb::core::bank::Sequence                  \n
     * Test of \ref gatb::core::bank::Sequence::getDataSize()   \n
     */
    void bank_checkMultipleFiles ()
    {
        /** We launch the test with uncompressed sample bank. */
        bank_checkMultipleFiles_aux (DBPATH("query.fa"));

        /** We launch the test with compressed sample bank. */
        bank_checkMultipleFiles_aux (DBPATH("query.fa.gz"));
    }

    /********************************************************************************/
    void bank_checkConvertBinary_aux (const string& filename, bool checkSize)
    {
        string filenameBin = filename + ".bin";

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare the two banks. */
        BankFasta       bank1 (filename);
        BankBinary bank2 (filenameBin);

        /** We convert the fasta bank in binary format. */
        BankFasta::Iterator itSeq1 (bank1);
        for (itSeq1.first(); !itSeq1.isDone(); itSeq1.next())   {  bank2.insert (*itSeq1);  }   bank2.flush ();

        /** We check that the binary bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filenameBin) == true);

        /** We check that the binary file size is smaller than the original file (should be about 25%). */
        if (checkSize)  {  CPPUNIT_ASSERT (System::file().getSize (filename) >= 4 * System::file().getSize (filenameBin));  }

        /** We suppress the binary file and check that the binary bank doesn't exist any more. */
        CPPUNIT_ASSERT (System::file().remove (filenameBin) == 0);
        CPPUNIT_ASSERT (System::file().doesExist (filenameBin) == false);
    }

    /** \brief Check that multiple files reading works well.
     */
    void bank_checkConvertBinary ()
    {
        bank_checkConvertBinary_aux (DBPATH("reads1.fa"),     true);
        bank_checkConvertBinary_aux (DBPATH("reads1.fa.gz"),  false);
        bank_checkConvertBinary_aux (DBPATH("reads2.fa"),     true);
    }

    /********************************************************************************/
    /** \brief Performance test
     */
    void bank_perf1 ()
    {
        const char* filename = "/media/GATB/users/edrezen/databases/swissprot.fa";
        //const char* filename = "/tmp/swissprot.fa.gz";
        //const char* filename = "/media/GATB/users/edrezen/databases/refseq_protein.00.fa";
        //const char* filename = "/media/GATB/users/edrezen/databases/uniprot.fa";

        /** We first clear the file system cache. */
        System::file().clearCache ();

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b (filename);

        /** We create an iterator over this bank. */
        BankFasta::Iterator it (b);

        ITime::Value t0 = System::time().getTimeStamp();

        u_int64_t totalSize = 0;

        /** We loop over sequences. */
        for (it.first(); !it.isDone(); it.next())
        {
            totalSize += it->getDataSize();
        }

        ITime::Value t1 = System::time().getTimeStamp();

        printf ("size=%lld  time=%lld => rate=%.3f MB/s \n",
            totalSize, t1-t0,
            (double)totalSize / (double) (t1-t0) / 1024 / 1024 * 1000
        );
    }

    /********************************************************************************/
    void bank_checkRegistery_aux (const string& bankformat, const string& bankuri, size_t nbCheck)
    {
        /** We get the default factory. */
        IBankFactory* factory = BankRegistery::singleton().getFactory (bankformat);

        /** We create a bank handle. */
        IBank* bank = factory->createBank (DBPATH(bankuri));
        LOCAL (bank);

        /** We create a bank iterator. */
        Iterator<Sequence>* itSeq = bank->iterator();
        LOCAL (itSeq);

        u_int64_t nbSeq = 0;

        /** We loop over sequences. */
        for (itSeq->first(); !itSeq->isDone(); itSeq->next())  {  nbSeq++;  }

        CPPUNIT_ASSERT (nbSeq == nbCheck);
    }

    /********************************************************************************/
    void bank_checkRegistery1 ()
    {
        bank_checkRegistery_aux (BankFasta::name(), "sample1.fa", 20);
    }

    /********************************************************************************/
    // We a define a functor that will be called during iteration for filtering some items.
    struct FilterFunctor  {  bool operator ()  (Sequence& seq)   {  return seq.getIndex() % 2 == 0; } };

    void bank_checkRegistery2 ()
    {
        /** We register our custom bank format. */
        BankRegistery::singleton().registerFactory (
            "customfasta",
            new BankFilteredFactory<FilterFunctor> ("fasta", FilterFunctor())
        );

        /** We should have only half the number of sequences of the original database. */
        bank_checkRegistery_aux ("customfasta", "sample1.fa", 10);
    }

    /********************************************************************************/
    void bank_strings1_aux (const char* table[], size_t actualNumber)
    {
        u_int64_t actualTotalSize=0, actualMaxSize=0;
        for (size_t i=0; i<actualNumber; i++)
        {
            actualTotalSize += strlen (table[i]);
            if (actualMaxSize < strlen(table[i]))  { actualMaxSize = strlen(table[i]); }
        }

        /** We create the bank. */
        BankStrings bank (table, actualNumber);

        Iterator<Sequence>* it = bank.iterator();
        LOCAL (it);

        size_t idx=0;
        for (it->first(); !it->isDone(); it->next(), idx++)
        {
            CPPUNIT_ASSERT ((*it)->getIndex()    == idx);
            CPPUNIT_ASSERT ((*it)->getDataSize() == strlen(table[idx]) );

            /** NOTE ! don't use strcmp since we consider buffers without ending '\0' */
            CPPUNIT_ASSERT (memcmp ((*it)->getDataBuffer(), table[idx], (*it)->getDataSize()) == 0);
        }

        CPPUNIT_ASSERT (idx == actualNumber);

        u_int64_t number, totalSize, maxSize;
        bank.estimate (number, totalSize, maxSize);
        CPPUNIT_ASSERT (number    == actualNumber);
        CPPUNIT_ASSERT (totalSize == actualTotalSize);
        CPPUNIT_ASSERT (maxSize   == actualMaxSize);
    }

    /********************************************************************************/
    void bank_strings1 ()
    {
        const char* table1[] =  {  "ACTACGATCGATGTA",  "TTAGAGCAGCGAG",  "AGGGGCCCATTTCATCTATC" };
        bank_strings1_aux (table1, ARRAY_SIZE (table1));

        const char* table2[] =  {  "ACTACGATCGATGTATTAGAGCAGCGAGAGGGGCCCATTTCATCTATC" };
        bank_strings1_aux (table2, ARRAY_SIZE (table2));

        const char* table3[] =  {
            "GATCCTCCCCAGGCCCCTACACCCAATGTGGAACCGGGGTCCCGAATGAAAATGCTGCTGTTCCCTGGAGGTGTTTTCCT",
            "GGACGCTCTGCTTTGTTACCAATGAGAAGGGCGCTGAATCCTCGAAAATCCTGACCCTTTTAATTCATGCTCCCTTACTC",
            "ACGAGAGATGATGATCGTTGATATTTCCCTGGACTGTGTGGGGTCTCAGAGACCACTATGGGGCACTCTCGTCAGGCTTC",
            "CGCGACCACGTTCCCTCATGTTTCCCTATTAACGAAGGGTGATGATAGTGCTAAGACGGTCCCTGTACGGTGTTGTTTCT",
            "GACAGACGTGTTTTGGGCCTTTTCGTTCCATTGCCGCCAGCAGTTTTGACAGGATTTCCCCAGGGAGCAAACTTTTCGAT"
        };
        bank_strings1_aux (table3, ARRAY_SIZE (table3));
    }

    /********************************************************************************/
    void bank_sequence ()
    {
        char* buf = (char*)"ACTACGATCGATGTA";

        Sequence s1 (buf);
        CPPUNIT_ASSERT (s1.getDataSize() == 15);
        CPPUNIT_ASSERT (strcmp (s1.getDataBuffer(), buf) == 0);
    }

    /********************************************************************************/
    void bank_random_aux (size_t nbSequences, size_t length)
    {
        BankRandom bank (nbSequences, length);

        Iterator<Sequence>* it = bank.iterator ();
        LOCAL (it);

        size_t foundSequences = 0;
        for (it->first(); !it->isDone(); it->next(), foundSequences++)
        {
            CPPUNIT_ASSERT ((*it)->getDataSize() == length);
        }

        CPPUNIT_ASSERT (foundSequences == nbSequences);
    }

    /** */
    void bank_random_1 ()
    {
        size_t nbTable[]  = {1, 5, 10, 100};
        size_t lenTable[] = {10, 100, 1000, 10*1000, 100*1000, 1000*1000};

        for (size_t i=0; i<ARRAY_SIZE(nbTable); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(lenTable); j++)
            {
                bank_random_aux (nbTable[i], lenTable[j]);
            }
        }
    }

    /********************************************************************************/
    void bank_splitter_1 ()
    {
        size_t   readSize = 20;
        u_int8_t coverage = 1;
        size_t   overlap  = 5;

        BankSplitter bank (
            new BankStrings ("ATCCTCCCCAGGCCCCTACACCCAATGTGGAACCGGGGTCCCGAATGAAAATGCTGCTGTTCCCTGGAGGTGTTCT", NULL),
            readSize, overlap, coverage
        );

        const char* check[] = {
          // ATCCTCCCCAGGCCCCTACACCCAATGTGGAACCGGGGTCCCGAATGAAAATGCTGCTGTTCCCTGGAGGTGTTCT
            "ATCCTCCCCAGGCCCCTACA",
                           "CTACACCCAATGTGGAACCG",
                                          "AACCGGGGTCCCGAATGAAA",
                                                         "TGAAAATGCTGCTGTTCCCT",
                                                                        "TCCCTGGAGGTGTTCT"
          // 0000000000111111111122222222223333333333444444444455555555556666666666777777
          // 0123456789012345678901234567890123456789012345678901234567890123456789012345
        };

        size_t nbIter = 0;
        Iterator<Sequence>* it = bank.iterator();  LOCAL (it);
        for (it->first(); !it->isDone(); it->next(), nbIter++)
        {
            string current = string ((*it)->getDataBuffer(), (*it)->getDataSize());
            CPPUNIT_ASSERT (current == check[nbIter]);
        }

        CPPUNIT_ASSERT (nbIter == ARRAY_SIZE(check));
    }

    /********************************************************************************/
    void bank_composite ()
    {
        const char* table[] =  {  "ACTACGATCGATGTA",  "TTAGAGCAGCGAG",  "AGGGGCCCATTTCATCTATC" };

        /** We create a composite bank. */
        vector<IBank*> banks;
        for (size_t i=0; i<ARRAY_SIZE(table); i++)  {  banks.push_back (new BankStrings (table[i], 0));  }

        /** We iterate the composite bank. */
        BankComposite bankComposite (banks);
        Iterator<Sequence>* it = bankComposite.iterator();
        LOCAL (it);

        size_t i=0;
        for (it->first(); !it->isDone(); it->next())
        {
            string data = string ((*it)->getDataBuffer(), (*it)->getDataSize());
            CPPUNIT_ASSERT ( data == table[i++]);
        }
        CPPUNIT_ASSERT (i==ARRAY_SIZE(table));
        CPPUNIT_ASSERT (bankComposite.getNbItems()==ARRAY_SIZE(table));
    }

    /********************************************************************************/
    void bank_album1 (void)
    {
        string albumUri = System::file().getTemporaryDirectory() + "/test_album.txt";

        /** We remove the album file (just to be sure). */
        System::file().remove (albumUri);

        /** We create an album. */
        BankAlbum album (albumUri);

        CPPUNIT_ASSERT (System::file().doesExist(albumUri) == true);

        /** We add a few files in the album. */
        album.add (System::file().getRealPath(DBPATH("reads1.fa")));
        album.add (System::file().getRealPath(DBPATH("reads2.fa")));

        /** We iterate the sequences. */
        Iterator<Sequence>* itSeq = album.iterator();
        LOCAL (itSeq);

        size_t nbSeq = 0;
        for (itSeq->first(); !itSeq->isDone(); itSeq->next())  { nbSeq++; }

        CPPUNIT_ASSERT (nbSeq == (100 + 1000) );  // reads1.fa has 100 seq, reads2.fa has 1000 seq

        /** We remove the album file. */
        System::file().remove (albumUri);

        CPPUNIT_ASSERT (System::file().doesExist(albumUri) == false);
    }

    /********************************************************************************/
    void bank_album2 (void)
    {
        /** We get the "album.txt" uri. */
        string albumUri = System::file().getRealPath(DBPATH("album.txt"));

        /** We check that the file exists. */
        CPPUNIT_ASSERT (System::file().doesExist(albumUri) == true);

        /** We create an album. */
        BankAlbum album (albumUri);

        /** We iterate the sequences. */
        Iterator<Sequence>* itSeq = album.iterator();
        LOCAL (itSeq);

        size_t nbSeq = 0;
        for (itSeq->first(); !itSeq->isDone(); itSeq->next())  { nbSeq++; }

        CPPUNIT_ASSERT (nbSeq == (100 + 1000) );  // reads1.fa has 100 seq, reads2.fa has 1000 seq
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestBank);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestBank);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

