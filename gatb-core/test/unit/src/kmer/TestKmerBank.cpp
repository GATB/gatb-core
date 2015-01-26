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
    template<class ModelType>
    void kmerbank_checkKmersFromBankAndBankBinary_aux (const char* filepath, size_t span)
    {
        /** Shortcuts. */
        typedef typename ModelType::Kmer     Kmer;
        typedef typename ModelType::Iterator ModelIterator;

        string filename    = DBPATH (filepath);
        string filenameBin = filepath + string(".bin");

        /** We check that the bank exists. */
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        /** We declare a kmer model with a given span size. */
        ModelType model (span);

        /** We declare the two banks.
         *  WARNING! we don't remove the bad characters (param 'false') in the binary bank,
         *  otherwise the two bank won't be comparable anymore. */
        BankFasta  bank1 (filename);
        BankBinary bank2 (filenameBin, false);

        /** We convert the fasta bank in binary format. */
        Iterator<Sequence>* itSeq1 = bank1.iterator();
        for (itSeq1->first(); !itSeq1->isDone(); itSeq1->next())   {  bank2.insert (itSeq1->item());  }   bank2.flush ();

        /** We declare two kmer iterators for the two banks and a paired one that links them. */
        ModelIterator* itKmer1 = new ModelIterator (model);   LOCAL (itKmer1);
        ModelIterator* itKmer2 = new ModelIterator (model);   LOCAL (itKmer2);

        PairedIterator<Kmer, Kmer> itKmer (itKmer1, itKmer2);

        /** We loop the two banks with a paired iterator. */
        Iterator<Sequence>* itSeq2 = bank2.iterator();
        PairedIterator<Sequence,Sequence> itSeq (itSeq1, itSeq2);

        /** We loop the sequences of the two banks. */
        for (itSeq.first(); !itSeq.isDone();  itSeq.next())
        {
            /** We set the data from which we want to extract kmers. */
            itKmer1->setData ((*itSeq1)->getData());
            itKmer2->setData ((*itSeq2)->getData());

            /** We loop the kmers for the two datas. */
            for (itKmer.first(); !itKmer.isDone();  itKmer.next())
            {
                CPPUNIT_ASSERT (itKmer->first.value() == itKmer->second.value());
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

        for (size_t i=0; i<ARRAY_SIZE(files); i++)
        {
            for (size_t j=0; j<ARRAY_SIZE(spans); j++)
            {
                kmerbank_checkKmersFromBankAndBankBinary_aux <Kmer<>::ModelDirect>    (files[i], spans[j]);
                kmerbank_checkKmersFromBankAndBankBinary_aux <Kmer<>::ModelCanonical> (files[i], spans[j]);
            }
        }
    }

    /********************************************************************************/
    void kmers_bankiterate ()
    {
        for (size_t i=0; i<12; i++)
        {
            u_int64_t nbKmers = 0;

            Iterator<Sequence>* it = BankKmers (i).iterator();  LOCAL (it);

            for (it->first(); !it->isDone(); it->next())  { nbKmers ++; }

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

