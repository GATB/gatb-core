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
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <iostream>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

#define DBPATH(a)  (string("../test/db/") + string(a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestKmerBank : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestKmerBank);

        CPPUNIT_TEST_GATB (kmerbank_checkKmersFromBankAndBankBinary);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    void kmerbank_checkKmersFromBankAndBankBinary_aux (const char* filepath, size_t span)
    {
        string filename    = DBPATH (filepath);
        string filenameBin = filepath + string(".bin");

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a kmer model with a given span size. */
        Model<u_int64_t> model (span);

        /** We declare the two banks. */
        Bank       bank1 (filename);
        BankBinary bank2 (filenameBin);

        /** We convert the fasta bank in binary format. */
        Bank::Iterator itSeq1 (bank1);
        for (itSeq1.first(); !itSeq1.isDone(); itSeq1.next())   {  bank2.insert (*itSeq1);  }   bank2.flush ();

        /** We declare two kmer iterators for the two banks and a paired one that links them. */
        Model<u_int64_t>::Iterator itKmer1 (model);
        Model<u_int64_t>::Iterator itKmer2 (model);
        PairedIterator<u_int64_t,u_int64_t> itKmer (itKmer1, itKmer2);

        /** We loop the two banks with a paired iterator. */
        BankBinary::Iterator itSeq2 (bank2);
        PairedIterator<Sequence,Sequence> itSeq (itSeq1, itSeq2);

        /** We loop the sequences of the two banks. */
        for (itSeq.first(); !itSeq.isDone();  itSeq.next())
        {
            /** We set the data from which we want to extract kmers. */
            itKmer1.setData (itSeq1->getData());
            itKmer2.setData (itSeq2->getData());

            /** We loop the kmers for the two datas. */
            for (itKmer.first(); !itKmer.isDone();  itKmer.next())
            {
                CPPUNIT_ASSERT (itKmer->first == itKmer->second);
            }
        }
    }

    /** \brief check that Bank and BankBinary provide the same set of kmers
     *
     * This test builds a binary bank from a fasta bank and check that the all the kmers from
     * initial bank are the same as the binary bank.
     */
    void kmerbank_checkKmersFromBankAndBankBinary ()
    {
        const char* files[] = { "reads1.fa", "reads1.fa.gz", "reads2.fa" };
        size_t      spans[] = { 2, 3, 5, 8, 13, 21, 34, 55 };

        for (size_t i=0; i<sizeof(files)/sizeof(files[0]); i++)
        {
            for (size_t j=0; j<sizeof(spans)/sizeof(spans[0]); j++)
            {
                kmerbank_checkKmersFromBankAndBankBinary_aux (files[i], spans[j]);
            }
        }
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION (TestKmerBank);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

