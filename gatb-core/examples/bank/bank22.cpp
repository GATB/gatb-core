//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*    Extract some sequences from a Fastq sequence file.                        */
/*                                                                              */
/* Cmd-line: bank22 <in_bank> <out_bank> <nb_seq_to_retain>                     */
/*                                                                              */
/* Sample: bank22 gatb-core/gatb-core/test/db/giab.hg002.2D.fastq.gz \          */
/*                /tmp/new_file.fastq                                           */
/*                500                                                           */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cerr << "you must provide: in_bank_name out_bank_name nb_seq_to_retain" << std::endl;
        return EXIT_FAILURE;
    }

	// We open the reference file
	IBank* inBank = Bank::open (argv[1]);
	std::string outB(argv[2]);
	//constructor: fileName to create,
	BankFasta outputBank (
		outB, // file name to create
		true  // write fastq instead of default fasta
		//,true // optional you can use gzip compression directly
	);

	// We create iterators over this bank.
	Iterator<Sequence>* itseq = inBank->iterator();
	itseq = inBank->iterator();
	LOCAL(itseq);
	int limit = atoi(argv[3]);
	int count=0;
	for (itseq->first(); !itseq->isDone(); itseq->next())
	{
		outputBank.insert (itseq->item());
		count++;
		if (count>=limit){
			break;
		}
	}
	outputBank.flush();
}
//! [snippet1]
