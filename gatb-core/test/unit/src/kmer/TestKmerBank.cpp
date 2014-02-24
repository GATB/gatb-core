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
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BankKmers.hpp>
#include <gatb/bank/impl/Banks.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <iostream>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

extern std::string DBPATH (const string& a);

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
        CPPUNIT_TEST_GATB (kmers_bankiterate);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    void kmerbank_checkKmersFromBankAndBankBinary_aux (const char* filepath, size_t span, KmerMode mode)
    {
        string filename    = DBPATH (filepath);
        string filenameBin = filepath + string(".bin");

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a kmer model with a given span size. */
        Kmer<>::Model model (span);

        /** We declare the two banks.
         *  WARNING! we don't remove the bad characters (param 'false') in the binary bank,
         *  otherwise the two bank won't be comparable anymore. */
        BankFasta       bank1 (filename);
        BankBinary bank2 (filenameBin, false);

        /** We convert the fasta bank in binary format. */
        BankFasta::Iterator itSeq1 (bank1);
        for (itSeq1.first(); !itSeq1.isDone(); itSeq1.next())   {  bank2.insert (*itSeq1);  }   bank2.flush ();

        /** We declare two kmer iterators for the two banks and a paired one that links them. */
        Kmer<>::Model::Iterator itKmer1 (model);
        Kmer<>::Model::Iterator itKmer2 (model);
        PairedIterator<Iterator, Kmer<>::Type,Kmer<>::Type> itKmer (itKmer1, itKmer2);

        /** We loop the two banks with a paired iterator. */
        BankBinary::Iterator itSeq2 (bank2);
        PairedIterator<Iterator, Sequence,Sequence> itSeq (itSeq1, itSeq2);

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
        size_t      spans[] = { 2, 3, 5, 8, 13, 21 };
        KmerMode    modes[] = { KMER_DIRECT, KMER_REVCOMP, KMER_MINIMUM};

        for (size_t i=0; i<ARRAY_SIZE(files); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(spans); j++)
            {
                for (size_t k=0; k<ARRAY_SIZE(modes); k++)
                {
                    kmerbank_checkKmersFromBankAndBankBinary_aux (files[i], spans[j], modes[k]);
                }
            }
        }
    }

    /********************************************************************************/
    struct Fct_kmers_bankiterate
    {
        Fct_kmers_bankiterate (u_int64_t& nbKmers) : nbKmers(nbKmers) {}
        u_int64_t& nbKmers;
        void operator() (Sequence& seq)  {  nbKmers ++;  }
    };

    void kmers_bankiterate ()
    {
        for (size_t i=0; i<12; i++)
        {
            u_int64_t nbKmers = 0;
#ifdef WITH_LAMBDA_EXPRESSIONS
            BankKmers (i).iterate ([&] (Sequence& seq)   {  nbKmers++; });
#else
            BankKmers (i).iterate (Fct_kmers_bankiterate(nbKmers));
#endif
            CPPUNIT_ASSERT (nbKmers == ((u_int64_t)1<<(2*i)));
        }
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestKmerBank);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestKmerBank);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

