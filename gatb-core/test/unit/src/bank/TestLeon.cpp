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

#include <gatb/tools/compression/Leon.hpp>

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

extern std::string DBPATH (const string& a);

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for Leon compress/decompressor
 */
class TestLeon : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestLeon);

    CPPUNIT_TEST_GATB(bank_checkLeon1);
    CPPUNIT_TEST_GATB(bank_checkLeon2);
    CPPUNIT_TEST_GATB(bank_checkLeon3);
    CPPUNIT_TEST_GATB(bank_checkLeon4);
    CPPUNIT_TEST_GATB(bank_checkLeon5);
    CPPUNIT_TEST_GATB(bank_checkLeon6);
	
	//removed some large files from distrib
   // CPPUNIT_TEST_GATB(bank_checkLeon7);
   // CPPUNIT_TEST_GATB(bank_checkLeon8);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/

    /**
     * Compare the content of two banks for strict equality: same number
     * of sequences, sequences in same order, etc.
     *
     * */
    void bank_compare_banks_equality(IBank* bank1, IBank* bank2){
		// We create iterators over this bank.
		Iterator<Sequence>* itFas = bank1->iterator();
		Iterator<Sequence>* itLeon = bank2->iterator();
		u_int64_t nbSeqFas = 0;
		u_int64_t nbSeqLeon = 0;
		{
			LOCAL(itFas);
			LOCAL(itLeon);
			// We do not use estimate() methods. Instead, we count
			// exact number of sequences in both banks
			for (itFas->first(); !itFas->isDone(); itFas->next()){nbSeqFas++;}
			for (itLeon->first(); !itLeon->isDone(); itLeon->next()){nbSeqLeon++;}
			CPPUNIT_ASSERT(nbSeqFas==nbSeqLeon);
		}
		itFas = bank1->iterator();
		itLeon = bank2->iterator();
		LOCAL(itFas);
		LOCAL(itLeon);
		// We create a PairedIterator to go through both banks simultaneously
		PairedIterator<Sequence,Sequence> it (itFas, itLeon);

		nbSeqFas =0;
		for (it.first(); !it.isDone(); it.next())
		{
			nbSeqFas++;
			// check sequence comment for equality
			CPPUNIT_ASSERT(it->first.getComment().compare(it->second.getComment())==0);
			// check sequence letters for equality
			CPPUNIT_ASSERT(it->first.toString().compare(it->second.toString())==0);
			// check sequence quality for equality
			CPPUNIT_ASSERT(it->first.getQuality().compare(it->second.getQuality())==0);
		}
		CPPUNIT_ASSERT(nbSeqFas==nbSeqLeon);
    }

    // http://stackoverflow.com/a/7026414
    /*******************************************************************************
     * From a fastq, generate leon compressed file, then compare that file to
     * the Fastq reference. Allow to check that Leon compression is still ok
     * over releases of GATB-Core.
     *
     * LOSSLESS version
     * */
    void bank_leon_compress_and_compare (const std::string& fastqFile, const std::string& leonFile)
    {
		// STEP 1: compress the Fastq reference file

		// we prepare the Leon command-line
    	std::vector<char*>       leon_args;
    	std::vector<std::string> data = {
    			"-",
				"-c",
				"-file", fastqFile,
				"-lossless", // <-- LOSSLESS
				"-verbose","0",
				"-kmer-size", "31",
				"-abundance", "1"
    	};
		for(std::vector<std::string>::iterator loop = data.begin(); loop != data.end(); ++loop){
			leon_args.push_back(&(*loop)[0]);
		}

		// we start Leon compressor
		Leon().run(leon_args.size(), &leon_args[0]);

		// STEP 2: compare reference and compressed version

		// we open the files in read mode
		
        IBank* fasBank = Bank::open (fastqFile); //BankFasta
		IBank* leonBank = Bank::open (leonFile); //BankLeon

		bank_compare_banks_equality(fasBank, leonBank);
		 
    }

    /**
     * Run Leon compress on a FastQ file.
     *
     * Parameter 'mode' is one of: 0, 1, 2 (for -noqual, -noheader
     * or -seq-only, respectively) or 3 (-reads 1000; used to test
     * parallel compression).
     */
    void run_leon_compressor(std::string& fastqFile, int mode){
		// we prepare the Leon command-line
    	std::vector<char*>       leon_args;
    	std::vector<std::string> data = {
    			"-",
				"-c",
				"-file", fastqFile,
				"-verbose","0",
				"-kmer-size", "31",
				"-abundance", "1",
    	}; // do NOT anymore modify this list (unless you also update
    	   // reference files). If you need to add other Leon args, see
    	   // switch() below.
    	switch (mode){
    	case 0:
    		data.push_back("-noqual");
    		break;
    	case 1:
    		data.push_back("-noheader");
    		break;
    	case 2:
    		data.push_back("-seq-only");
    		break;
    	case 3:
    		data.push_back("-lossless");
    		data.push_back("-reads");
    		data.push_back("1000");
    		break;
    	}

		for(std::vector<std::string>::iterator loop = data.begin(); loop != data.end(); ++loop){
			leon_args.push_back(&(*loop)[0]);
		}

		// we start Leon compressor
		Leon().run(leon_args.size(), &leon_args[0]);

    }

    /**
     * Check the content of a Leon compressed file according to the use
     * of Leon's argument -noqual, -noheader or -seq-only (mode= 0, 1 or 2,
     * respectively).
     */
    void check_leon_content(std::string& leonFile, int mode){
    	// we open the leon file in read mode
		IBank* leonBank = Bank::open (leonFile); //BankLeon

		// We create iterators over this bank.
		Iterator<Sequence>* itLeon = leonBank->iterator();
		itLeon = leonBank->iterator();
		LOCAL(itLeon);
		for (itLeon->first(); !itLeon->isDone(); itLeon->next())
		{
			Sequence& seq = itLeon->item();
	    	switch (mode){
	    	case 0://"-noqual"
	    		CPPUNIT_ASSERT(seq.getComment().size()!=0);
				CPPUNIT_ASSERT(seq.toString().size()!=0);
				CPPUNIT_ASSERT(seq.getQuality().size()==0);
	    		break;
	    	case 1://"-noheader"
	    		// noheader: comment contains read rank number starting from 0
				CPPUNIT_ASSERT(std::stoi(seq.getComment())>=0);
				CPPUNIT_ASSERT(seq.toString().size()!=0);
				CPPUNIT_ASSERT(seq.getQuality().size()!=0);
	    		break;
	    	case 2://"-seq-only"
	    		// seq-only: comment contains read rank number starting from 0
				CPPUNIT_ASSERT(std::stoi(seq.getComment())>=0);
				CPPUNIT_ASSERT(seq.toString().size()!=0);
				CPPUNIT_ASSERT(seq.getQuality().size()==0);
	    		break;
	    	}
		}
    }

    /********************************************************************************
     * We compare a lossless-compressed Leon file against the original Fastq file.
     * Allow to ensure that Leon decompression is still working over GATB-Core
     * releases.
     */
    void bank_checkLeon1 ()
    {
    	// We open and test the reference file
    	string fasPath=DBPATH("leon1.fastq");
    	string btype = Bank::getType(fasPath);
    	CPPUNIT_ASSERT(
				btype.compare("fasta")==0 ||
				btype.compare("fastq")==0
		);
		IBank* fasBank = Bank::open (fasPath);

    	// We open and test the leon file
		// Caution: this must be a LOSSLESS-compressed file
		// (because we compare quality, see below)
    	string leonPath=DBPATH("leon1.fastq.leon-ref");
		btype = Bank::getType(leonPath);
		CPPUNIT_ASSERT(btype.compare("leon")==0);
		IBank* leonBank = Bank::open (leonPath);

		bank_compare_banks_equality(fasBank, leonBank);
    }

    /*******************************************************************************
     * From a fastq, generate leon compressed file, then compare that file to
     * the Fastq reference. Allow to check that Leon compression is still ok
     * over releases of GATB-Core.
     *
     * LOSSLESS version
     * */
    void bank_checkLeon2(){
    	bank_leon_compress_and_compare(DBPATH("leon2.fastq"), DBPATH("leon2.fastq.leon"));
    }

    /*******************************************************************************
     * From a fastq, generate leon compressed file, then compare that file to
     * the leon compressed reference. Allow to check that still format is still ok
     * over releases of GATB-Core.
     *
     * LOSSLESS version
     * */
    void bank_checkLeon3 ()
    {
    	// The existing reference file
    	std::string fastqFile = DBPATH("leon2.fastq");
    	// The Leon file to create
		string leonFile=fastqFile+".leon";
    	// The Leon file to use as a reference
		string leonFileRef=fastqFile+".leon-ref";

		// STEP 1: compress the Fasta file

		// we prepare the Leon command-line
    	std::vector<char*>       leon_args;
    	std::vector<std::string> data = {
    			"-",
				"-c",
				"-file", fastqFile,
				"-lossless", // <-- LOSSLESS
				"-verbose","0",
				"-kmer-size", "31",
				"-abundance", "1"
    	};
		for(std::vector<std::string>::iterator loop = data.begin(); loop != data.end(); ++loop){
			leon_args.push_back(&(*loop)[0]);
		}

		// we start Leon compressor
		Leon().run(leon_args.size(), &leon_args[0]);

		// STEP 2: compare reference and compressed version

		// we open the files in read mode
        IBank* leonBank = Bank::open (leonFile); //BankLeon freshly created
		IBank* leonBankRef = Bank::open (leonFileRef); //BankLeon reference

		bank_compare_banks_equality(leonBank, leonBankRef);
    }

    /*******************************************************************************
	 * Test other args of Leon:
	 *   -noqual
	 *
	 * */
	void bank_checkLeon4 ()
	{
    	// The existing reference file
    	std::string fastqFile = DBPATH("leon2.fastq");

    	// The Leon file to create
		string leonFile=fastqFile+".leon";

		// STEP 1: compress the Fasta file
		run_leon_compressor(fastqFile, 0);

		// STEP 2: check leon content
		check_leon_content(leonFile, 0);

	}

    /*******************************************************************************
	 * Test other args of Leon:
	 *   -noheader
	 *
	 * */
	void bank_checkLeon5 ()
	{
    	// The existing reference file
    	std::string fastqFile = DBPATH("leon2.fastq");

    	// The Leon file to create
		string leonFile=fastqFile+".leon";

		// STEP 1: compress the Fasta file
		run_leon_compressor(fastqFile, 1);

		// STEP 2: check leon content
		check_leon_content(leonFile, 1);

	}
    /*******************************************************************************
	 * Test other args of Leon:
	 *   -seq-only
	 *
	 * */
	void bank_checkLeon6 ()
	{
    	// The existing reference file
    	std::string fastqFile = DBPATH("leon2.fastq");

    	// The Leon file to create
		string leonFile=fastqFile+".leon";

		// STEP 1: compress the Fastq file
		run_leon_compressor(fastqFile, 2);

		// STEP 2: check leon content
		check_leon_content(leonFile, 2);
	}

	/**
	 * Same as bank_checkLeon2() but with a bigger file.
	 * */
	void bank_checkLeon7 ()
	{
		bank_leon_compress_and_compare(
				DBPATH("NIST7035_TAAGGCGA_L001_R1_001_5OK.fastq.gz"),
				DBPATH("NIST7035_TAAGGCGA_L001_R1_001_5OK.fastq.leon")
				);
	}

    /*******************************************************************************
	 * Test Leon parallel compression/decompression.
	 *
	 * */
	void bank_checkLeon8 ()
	{
    	// The existing FastQ reference file (contains 50K reads)
		std::string fileName =  DBPATH("NIST7035_TAAGGCGA_L001_R1_001_5OK.fastq");
		std::string fastqFile = fileName + ".gz";
    	// The Leon file to create (WILL BE compressed NOW using -reads 1000)
		string leonFile = fileName + ".leon";
    	// The Leon reference file (WAS compressed using -reads 50000)
		string leonRefFile = leonFile+"-ref";

		// STEP 1: compress the Fastq file ('mode=3' means use '-reads 1000')
		run_leon_compressor(fastqFile, 3);

		// STEP 2: compare reference and compressed version
        IBank* leonRefBank = Bank::open (leonRefFile);
		IBank* leonBank = Bank::open (leonFile);
		bank_compare_banks_equality(leonRefBank, leonBank);
	}
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestLeon);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestLeon);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

