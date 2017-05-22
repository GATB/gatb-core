//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/* !!!!!  WARNING !!!!!
 * DO NOT EDIT: snippet used to test Leon!
 */
#undef NDEBUG
#include <assert.h>

/********************************************************************************/
/*                    Compare two banks for equality                            */
/*                                                                              */
/* This snippet shows how to compare two banks:                                 */
/*   - first bank is either a Fasta or a Fastq file                             */
/*   - secon bank is either a Fasta or a Fastq or a Leon lossless file          */
/*                                                                              */
/* Note: we use here a try/catch block in case the bank opening doesn't work.   */
/*                                                                              */
/* Cmd-line: bank25 <fasta/q file> <leon file>                                  */
/*                                                                              */
/* Sample: bank25 gatb-core/gatb-core/test/db/sample.fastq \                    */
/*                gatb-core/gatb-core/test/db/sample.fastq.leon                 */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "you must provide two banks." << std::endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        string btype = Bank::getType(argv[1]);
    	assert(
       		btype.compare("fasta")==0 ||
		btype.compare("fastq")==0
        );
        // We open the reference file
        IBank* fasBank = Bank::open (argv[1]);

        btype = Bank::getType(argv[2]);

        // We open its leon-lossless-compressed representation
        IBank* leonBank = Bank::open (argv[2]);

        u_int64_t nbSeqFas = 0;
        u_int64_t nbSeqLeon = 0;

        // We create iterators over this bank.
        Iterator<Sequence>* itFas = fasBank->iterator();
        Iterator<Sequence>* itLeon = leonBank->iterator();
        {
		// we use a GATB-Core macro to automatically release
		// memory allocated here as soon as we leave this code
		// block
		LOCAL(itFas);
		LOCAL(itLeon);

		// We do not use estimate() methods. Instead, we count
		// exact number of sequences in both banks
		for (itFas->first(); !itFas->isDone(); itFas->next()){nbSeqFas++;}
		for (itLeon->first(); !itLeon->isDone(); itLeon->next()){nbSeqLeon++;}
		assert(nbSeqFas==nbSeqLeon);
        }

	// We create a PairedIterator to go through both banks simultaneously
	itFas = fasBank->iterator();
	itLeon = leonBank->iterator();
	LOCAL(itFas);
	LOCAL(itLeon);
	PairedIterator<Sequence,Sequence> it (itFas, itLeon);

        for (it.first(); !it.isDone(); it.next())
        {
        	// check sequence comment for equality
        	assert(it->first.getComment().compare(it->second.getComment())==0);
        	// check sequence letters for equality
        	assert(it->first.toString().compare(it->second.toString())==0);
        	// check sequence quality for equality
        	assert(it->first.getQuality().compare(it->second.getQuality())==0);
        }
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
