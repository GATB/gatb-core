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

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankStrings.hpp>

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
        //CPPUNIT_TEST_GATB (bank_perf1);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/
    void bank_checkSample1_aux (const string& filename, Bank::Iterator::CommentMode_e mode)
    {
        const char* text = "ARNDCQEGHILKMFPSTWYV";

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        Bank b (filename);

        /** We create an iterator over this bank. */
        Bank::Iterator it (b, mode);

        size_t i=0;

        /** We loop over sequences. */
        for (it.first(); !it.isDone(); it.next(), i++)
        {
            /** We prepare the sequence string to be matched according to the provided mode. */
            char buffer[32];
            snprintf (buffer, sizeof(buffer), "seq%d%s", i+1,  (mode == Bank::Iterator::FULL ? " generic" : ""));

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
        bank_checkSample1_aux (DBPATH("sample1.fa"), Bank::Iterator::IDONLY);
        bank_checkSample1_aux (DBPATH("sample1.fa"), Bank::Iterator::FULL);

        /** We launch the test with compressed sample bank. */
        bank_checkSample1_aux (DBPATH("sample1.fa.gz"), Bank::Iterator::IDONLY);
        bank_checkSample1_aux (DBPATH("sample1.fa.gz"), Bank::Iterator::FULL);
    }

    /********************************************************************************/
    void bank_checkSample2_aux (const string& filename)
    {
        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        Bank b (filename);

        /** We create an iterator over this bank with only id comments. */
        Bank::Iterator it (b, Bank::Iterator::IDONLY);

        size_t nbSeq=0;

        /** We loop over sequences. */
        for (it.first(); !it.isDone(); it.next(), nbSeq++)
        {
            char buffer[32];
            snprintf (buffer, sizeof(buffer), "seq%d", nbSeq+1);

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
        Bank b1 (filename);

        /** We iterate without comments. */
        for (Bank::Iterator it (b1, Bank::Iterator::NONE); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty());
        }

        /** We iterate with only id comments. */
        for (Bank::Iterator it (b1, Bank::Iterator::IDONLY); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty() == false);
            CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") == 0);
        }

        /** We iterate with only full comments. */
        for (Bank::Iterator it (b1, Bank::Iterator::FULL); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->getComment().empty() == false);
            CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") != 0);
        }

        /** We iterate with default value (should be FULL). */
        for (Bank::Iterator it (b1); !it.isDone(); it.next())
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
        Bank bankKO ("_dummy_bank_name_");
        Bank bankOK (DBPATH("sample1.fa"));

        /** We check that we can't get a valid iterator on it. */
        CPPUNIT_ASSERT_THROW (Bank::Iterator it(bankKO), gatb::core::system::Exception);

        /** We check that the bank does exist. */
        CPPUNIT_ASSERT_NO_THROW (Bank::Iterator it(bankOK));
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
        Bank b1 (DBPATH("sample1.fa"));

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
        CPPUNIT_ASSERT_THROW (Bank b (filenames), gatb::core::system::Exception);

        while (filenames.size() < Bank::getMaxNbFiles())
        {
            filenames.push_back(DBPATH("sample1.fa"));

            CPPUNIT_ASSERT_NO_THROW (Bank b (filenames));
        }

        /** We check that we can't use a Bank with with too many filenames. */
        filenames.push_back(DBPATH("sample1.fa"));
        CPPUNIT_ASSERT_THROW (Bank b (filenames), gatb::core::system::Exception);
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
        Bank b1 (DBPATH("sample1.fa"));

        /** We check the estimation of sequences number. */
        u_int64_t estim1 = b1.estimateNbSequences();
        CPPUNIT_ASSERT (estim1 == 20);

        /** We build another bank holding several time the same bank. */
        for (size_t i=1; i<=Bank::getMaxNbFiles(); i++)
        {
            vector<string> filenames;
            for (size_t j=1; j<=i; j++)  { filenames.push_back(DBPATH("sample1.fa")); }

            Bank b2 (filenames);

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
    void bank_checkProgress_aux (Bank& b, size_t modulo)
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
        Bank b1 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b1, tableMod[i]);  }

        /** We build a bank holding x2 the same bank.*/
        while (filenames.size() < 2)   {  filenames.push_back (filename);  }
        Bank b2 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b2, tableMod[i]);  }

        /** We build a bank holding x10 the same bank.*/
        while (filenames.size() < 10)   {  filenames.push_back (filename);  }
        Bank b10 (filenames);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b10, tableMod[i]);  }

        /** We build a bank holding x20 the same bank.*/
        while (filenames.size() < 20)   {  filenames.push_back (filename);  }
        Bank b20 (filenames);
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
        Bank b1 (filename);

        /** We gather some information about this (single) bank. */
        Info i1;
        for (Bank::Iterator it (b1); !it.isDone(); it.next())   {  i1.nbseq++;   i1.datasize += it->getDataSize(); }

        /** We build a bank holding x2 the same bank.*/
        while (filenames.size() < 2)   {  filenames.push_back (filename);  }
        Bank b2 (filenames);

        /** We gather some information about this x2 bank. */
        Info i2;
        for (Bank::Iterator it (b2); !it.isDone(); it.next())   {  i2.nbseq++;   i2.datasize += it->getDataSize(); }

        /** We build a bank holding x4 the same bank.*/
        while (filenames.size() < 4)   {  filenames.push_back (filename);  }
        Bank b4 (filenames);

        /** We gather some information about this x4 bank. */
        Info i4;
        for (Bank::Iterator it (b4); !it.isDone(); it.next())   {  i4.nbseq++;   i4.datasize += it->getDataSize(); }

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
        Bank       bank1 (filename);
        BankBinary bank2 (filenameBin);

        /** We convert the fasta bank in binary format. */
        Bank::Iterator itSeq1 (bank1);
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
        Bank b (filename);

        /** We create an iterator over this bank. */
        Bank::Iterator it (b);

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
        bank_checkRegistery_aux (Bank::name(), "sample1.fa", 20);
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
            CPPUNIT_ASSERT (strcmp ((*it)->getDataBuffer(), table[idx]) == 0);
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

};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION (TestBank);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

