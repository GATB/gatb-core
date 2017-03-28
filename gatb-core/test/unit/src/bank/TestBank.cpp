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

#include <gatb/bank/impl/Bank.hpp>
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

using namespace gatb::core::tools::misc;

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
        CPPUNIT_TEST_GATB (bank_checkSample3);
        CPPUNIT_TEST_GATB (bank_checkComments);
        CPPUNIT_TEST_GATB (bank_checkBadUri);
        CPPUNIT_TEST_GATB (bank_checkSize);
        CPPUNIT_TEST_GATB (bank_checkMultipleFiles);
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
        CPPUNIT_TEST_GATB (bank_album3);
        CPPUNIT_TEST_GATB (bank_iteration);
        //        CPPUNIT_TEST_GATB (bank_datalinesize); // disabled since we're printing fasta in one line now (see "#if 1" in BankFasta)
        CPPUNIT_TEST_GATB (bank_registery_types);
        CPPUNIT_TEST_GATB (bank_checkPower2);

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
    void bank_checkSample3_aux (const string& filename)
    {
        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        IBank* bank = Bank::open(filename);
        CPPUNIT_ASSERT (bank != NULL);
        LOCAL (bank);

        /** We create an iterator over this bank. */
        Iterator<Sequence>* it = bank->iterator();
        CPPUNIT_ASSERT (it != NULL);
        LOCAL (it);

        size_t nbSeq=0;

        /** We loop over sequences. */
        for (it->first(); !it->isDone(); it->next(), nbSeq++)
        {
            CPPUNIT_ASSERT (it->item().getDataSize() > 0);
        }

        /** We check that we read exactly N sequences. */
        CPPUNIT_ASSERT (nbSeq == 7);
    }
    /********************************************************************************/
    void bank_checkSample3 ()
    {
        /** We launch the test with fastq banks. */
        bank_checkSample3_aux (DBPATH("sample.fastq"));
        bank_checkSample3_aux (DBPATH("sample.fastq.gz"));
    }

    /********************************************************************************/
    void bank_checkComments_aux (const string& filename)
    {
        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b1 (filename);

        /** We iterate without comments. */
        {
            BankFasta::Iterator it (b1, BankFasta::Iterator::NONE);
            for (it.first(); !it.isDone(); it.next())
            {
                CPPUNIT_ASSERT (it->getComment().empty());
            }
        }

        /** We iterate with only id comments. */
        {
            BankFasta::Iterator it (b1, BankFasta::Iterator::IDONLY);
            for (it.first(); !it.isDone(); it.next())
            {
                CPPUNIT_ASSERT (it->getComment().empty() == false);
                CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") == 0);
            }
        }

        /** We iterate with only full comments. */
        {
            BankFasta::Iterator it (b1, BankFasta::Iterator::FULL);
            for (it.first(); !it.isDone(); it.next())
            {
                CPPUNIT_ASSERT (it->getComment().empty() == false);
                CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") != 0);
            }
        }

        /** We iterate with default value (should be FULL). */
        {
            BankFasta::Iterator it (b1);
            for (it.first(); !it.isDone(); it.next())
            {
                CPPUNIT_ASSERT (it->getComment().empty()==false);
                CPPUNIT_ASSERT (strstr (it->getComment().c_str(), " ") != 0);
            }
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
    /** \brief Check that we can estimate the number of sequences in a bank
     *
     * Since we know the number of sequences in one sample bank, we can correctly
     * retrieve this number through the IBank API.
     *
     * The test is done with one simple sample, then with banks defined by several files.
     *
     * Test of \ref gatb::core::bank::impl::Bank                        \n
     * Test of \ref gatb::core::bank::impl::Bank::estimateNbSequences() \n
     */
    void bank_checkEstimateNbSequences ()
    {
        string filename = DBPATH("sample1.fa");

        /** We declare a Bank instance. */
        BankFasta b1 (filename);

        /** We check the estimation of sequences number. */
        u_int64_t estim1 = b1.estimateNbItems();
        CPPUNIT_ASSERT (estim1 == 20);

        string albumName = "album.txt";
        if (System::file().doesExist (albumName) == true)  { System::file().remove (albumName); }

        /** We create the album bank. */
        BankAlbum album (albumName);

        size_t nbMaxFiles = 30;

        /** We build another bank holding several time the same bank. */
        for (size_t i=1; i<=nbMaxFiles; i++)
        {
            /** Add a bank to the album. */
            album.addBank (filename);

            CPPUNIT_ASSERT ((size_t)album.estimateNbItems() == i * estim1);
        }

        CPPUNIT_ASSERT (album.getNbBanks() == nbMaxFiles);

        CPPUNIT_ASSERT (System::file().remove (albumName) == 0);
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
    void bank_checkProgress_aux (IBank& b, size_t modulo)
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
        CPPUNIT_ASSERT (fct->nbCalls == (b.estimateNbItems() + modulo - 1) / modulo);

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
        BankFasta b1 (filename);
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (b1, tableMod[i]);  }

        /** We create the album bank. */
        string albumName = "album.txt";
        if (System::file().doesExist (albumName) == true)  { System::file().remove (albumName); }
        BankAlbum album (albumName);

        /** We build a bank holding x2 the same bank.*/
        while (album.getNbBanks() < 2)   {  album.addBank (filename);  }
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (album, tableMod[i]);  }

        /** We build a bank holding x10 the same bank.*/
        while (album.getNbBanks() < 10)   {  album.addBank (filename);  }
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (album, tableMod[i]);  }

        /** We build a bank holding x20 the same bank.*/
        while (album.getNbBanks() < 20)   {  album.addBank (filename);  }
        for (size_t i=0; i<sizeof(tableMod)/sizeof(tableMod[0]); i++)  {  bank_checkProgress_aux (album, tableMod[i]);  }

        System::file().remove (albumName);
    }

    /********************************************************************************/
    void bank_checkMultipleFiles_aux (const string& filename)
    {
        /** A utility structure for this test. */
        struct Info  {  u_int32_t nbseq;  u_int64_t datasize;    Info() : nbseq(0), datasize(0) {}  };

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a Bank instance. */
        BankFasta b1 (filename);

        /** We gather some information about this (single) bank. */
        Info i1;
        BankFasta::Iterator it (b1);
        for (it.first(); !it.isDone(); it.next())   {  i1.nbseq++;   i1.datasize += it->getDataSize(); }

        /** We create the album bank. */
        string albumName = "album.txt";
        if (System::file().doesExist (albumName) == true)  { System::file().remove (albumName); }
        BankAlbum album (albumName);

        /** We build a bank holding x2 the same bank.*/
        while (album.getNbBanks() < 2)   {  album.addBank (filename);  }
        CPPUNIT_ASSERT (album.getNbBanks() == 2);

        /** We gather some information about this x2 bank. */
        Info i2;
        Iterator<Sequence>* it2 = album.iterator();  LOCAL (it2);
        for (it2->first(); !it2->isDone(); it2->next())   {  i2.nbseq++;   i2.datasize += (*it2)->getDataSize(); }

        /** We build a bank holding x4 the same bank.*/
        while (album.getNbBanks() < 4)   {  album.addBank (filename);  }
        CPPUNIT_ASSERT (album.getNbBanks() == 4);

        /** We gather some information about this x4 bank. */
        Info i4;
        Iterator<Sequence>* it4 = album.iterator();  LOCAL (it4);
        for (it4->first(); !it4->isDone(); it4->next())   {  i4.nbseq++;   i4.datasize += (*it4)->getDataSize(); }

        CPPUNIT_ASSERT (i1.nbseq   * 2 == i2.nbseq);
        CPPUNIT_ASSERT (i1.datasize* 2 == i2.datasize);

        CPPUNIT_ASSERT (i1.nbseq   * 4 == i4.nbseq);
        CPPUNIT_ASSERT (i1.datasize* 4 == i4.datasize);

        CPPUNIT_ASSERT (i2.nbseq   * 2 == i4.nbseq);
        CPPUNIT_ASSERT (i2.datasize* 2 == i4.datasize);

        System::file().remove (albumName);
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
        BankFasta  bank1 (filename);
        BankBinary bank2 (filenameBin);

        /** We convert the fasta bank in binary format. */
        BankFasta::Iterator itSeq1 (bank1);
        for (itSeq1.first(); !itSeq1.isDone(); itSeq1.next())   {  bank2.insert (*itSeq1);  }   bank2.flush ();

        /** We check that the binary bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filenameBin) == true);

        /** We check that the binary file size is smaller than the original file (should be about 25%). */
        if (checkSize)  {  CPPUNIT_ASSERT (System::file().getSize (filename) >= 4 * System::file().getSize (filenameBin));  }

        Iterator<Sequence>* it1 = bank1.iterator();  LOCAL (it1);
        Iterator<Sequence>* it2 = bank2.iterator();  LOCAL (it2);

        /** We check that both banks have the same nucleotides. */
        PairedIterator<Sequence> itPair (it1, it2);
        for (itPair.first(); !itPair.isDone(); itPair.next())
        {
            Sequence& s1 = itPair.item().first;
            Sequence& s2 = itPair.item().second;

            CPPUNIT_ASSERT (s1.getDataSize() == s2.getDataSize());

            size_t len = s1.getDataSize();

            char* b1 = s1.getDataBuffer();
            char* b2 = s2.getDataBuffer();

            for (size_t i=0; i<len; i++)
            {
                Data::ConvertChar c1 = Data::ConvertASCII::get  (b1, i);
                Data::ConvertChar c2 = Data::ConvertBinary::get (b2, i);
                CPPUNIT_ASSERT (c1.first == c2.first);
            }
        }

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

        printf ("size=%ld  time=%ld => rate=%.3f MB/s \n",
            totalSize, t1-t0,
            (double)totalSize / (double) (t1-t0) / 1024 / 1024 * 1000
        );
    }

    /********************************************************************************/
    void bank_checkRegistery_aux (const string& bankformat, const string& bankuri, size_t nbCheck)
    {
        /** We create a bank handle. */
        IBank* bank = Bank::open (DBPATH(bankuri));
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
        string factoryName = "customfasta";

        /** We register our custom bank format. */
        Bank::registerFactory (
                factoryName,
            new BankFilteredFactory<FilterFunctor> ("fasta", FilterFunctor()),
            true  // => put this factory in first position to override fasta factory
        );

        /** We should have only half the number of sequences of the original database. */
        bank_checkRegistery_aux (factoryName, "sample1.fa", 10);

        /** We have to unregister our custom factory now (default fasta factory will be used again). */
        CPPUNIT_ASSERT (Bank::unregisterFactory (factoryName) == true);
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
        /** The idea is to create an album bank from scratch and add some banks into. */

        string albumUri = System::file().getTemporaryDirectory() + "/test_album.txt";

        /** We remove the album file (just to be sure). */
        System::file().remove (albumUri);

        /** We create an album. */
        BankAlbum album (albumUri);

        CPPUNIT_ASSERT (System::file().doesExist(albumUri) == true);

        /** We add a few files in the album. */
        album.addBank (System::file().getRealPath(DBPATH("reads1.fa")));
        album.addBank (System::file().getRealPath(DBPATH("reads2.fa")));

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

    /********************************************************************************/
    void bank_album3_aux (const string& filename, size_t nbBanks, size_t nbSeqCheck)
    {
        string albumFilename = "test.album";

        ofstream albumFile (albumFilename.c_str());
        CPPUNIT_ASSERT (albumFile.is_open());

        for (size_t i=0; i<nbBanks; i++)  {  albumFile << filename << endl; }

        IBank* bank = Bank::open (albumFilename);
        LOCAL (bank);

        /** We iterate the sequences. */
        Iterator<Sequence>* itSeq = bank->iterator();
        LOCAL (itSeq);

        size_t nbSeq = 0;
        for (itSeq->first(); !itSeq->isDone(); itSeq->next())  { nbSeq++; }

        CPPUNIT_ASSERT (nbSeq == nbSeqCheck*nbBanks);

        /** Cleanup. */
        System::file().remove(albumFilename);
    }

    void bank_album3 ()
    {
        bank_album3_aux (DBPATH("reads1.fa"),       2000, 100);
        bank_album3_aux (DBPATH("sample.fastq.gz"), 2000, 7);
        bank_album3_aux (DBPATH("sample.fastq"),    2000, 7);
    }

    /********************************************************************************/

    struct Fct { void operator() (Sequence& s) { cout << s.getComment() << endl;} };

    /********************************************************************************/
    void bank_iteration (void)
    {
        BankFasta bank (DBPATH("reads1.fa"));

        size_t count = 0;

        Iterator<Sequence>* it = bank.iterator();  LOCAL (it);
        for (it->first(); !it->isDone(); it->next())  { count ++; }
        CPPUNIT_ASSERT (count == 100);
    }

    /********************************************************************************/
    void bank_datalinesize_aux (const char* sequence, size_t dataLineSize)
    {
        /** We create a fake bank. */
        BankStrings inputBank (sequence, NULL);

        string outputFilename = System::file().getTemporaryDirectory() + "/foo.fa";
        System::file().remove (outputFilename);
        CPPUNIT_ASSERT (System::file().doesExist(outputFilename) == false);

        BankFasta::setDataLineSize (dataLineSize);

        /** We want the output bank instance to be deleted before reading the generated file. */
        {
            /** We create a FASTA bank. */
            BankFasta outputBank (outputFilename);

            /** We iterate the sequences. */
            Iterator<Sequence>* itSeq = inputBank.iterator();
            LOCAL (itSeq);

            for (itSeq->first(); !itSeq->isDone(); itSeq->next())  {  outputBank.insert (itSeq->item());  }
            outputBank.flush();
        }

        const char* loop = sequence;

        ifstream fbank (outputFilename.c_str());
        if (fbank.is_open())
        {
            string line;
            while (getline (fbank,line))
            {
                /** We skip the comment. */
                if (line[0] == '>')  { continue; }

                CPPUNIT_ASSERT (line.size() == dataLineSize);

                for (size_t i=0; i<dataLineSize; i++, loop++)
                {
                    CPPUNIT_ASSERT (line[i] == *loop);
                }
            }
            fbank.close();
        }
        else
        {
            CPPUNIT_ASSERT (false);
        }

        /** We remove the album file. */
        System::file().remove (outputFilename);
        CPPUNIT_ASSERT (System::file().doesExist(outputFilename) == false);
    }

    /********************************************************************************/
    void bank_datalinesize (void)
    {
        // Important : sequence size power of 2 for the test
        const char* sequence = "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACCGGGTATATGGCGCGCGATTATATACGCGCGCTATATACGCGCTATCGATCGATCGAGCGACTAAT";

        size_t dataLineSizeTable[] = {1, 2, 4, 8, 16, 32, 64, 128};

        for (size_t i=0; i<ARRAY_SIZE(dataLineSizeTable); i++)
        {
            bank_datalinesize_aux (sequence, dataLineSizeTable[i]);
        }
    }

    /********************************************************************************/
    void bank_registery_types (void)
    {
        IBank* bank1 = Bank::open (DBPATH("sample1.fa"));
        CPPUNIT_ASSERT (bank1 != 0);
        LOCAL (bank1);

        CPPUNIT_ASSERT (Bank::getType(DBPATH("album.txt"))    == "album");
        CPPUNIT_ASSERT (Bank::getType(DBPATH("sample1.fa"))   == "fasta");
        CPPUNIT_ASSERT (Bank::getType(DBPATH("sample.fastq")) == "fastq");
    }

    /********************************************************************************/
    void bank_checkPower2 ()
    {
        string filename = "test.fa";

        /** We create a file with a specific content (ie a space at idx being a power of 2, which used make
         * and invalid write detected by valgrind). */
        ofstream file (filename.c_str());
        CPPUNIT_ASSERT (file.is_open());

        file << ">123456789012345 789012345678901" << endl;
        file << "CTTGTAATAACATGCACAAATACAATGCCATTCAATTTGA";
        file.flush ();

        /** We open it as a fasta bank. */
        BankFasta b (filename);

        BankFasta::Iterator it (b);
        for (it.first(); !it.isDone(); it.next())  {  CPPUNIT_ASSERT (it->getDataSize() == 40);  }

        /** Some cleanup. */
        System::file().remove(filename);
        CPPUNIT_ASSERT (System::file().doesExist(filename) == false);
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestBank);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestBank);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

