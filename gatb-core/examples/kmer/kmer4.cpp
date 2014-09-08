//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <string.h>

// We use the required packages
using namespace std;

static const size_t span = KSIZE_1;

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how to iterate the kmers of the sequences from a bank.    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"           << endl;
        cerr << "   2) list of FASTA files" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        u_int64_t nbSequences = 0;
        u_int64_t nbKmers     = 0;

        // We declare a Bank instance defined by a list of filenames
        BankFasta b (argc-2, argv+2);

        // We declare a kmer model with a given span size.
        Kmer<span>::ModelDirect model (kmerSize);

        // We create an iterator over this bank.
        BankFasta::Iterator itSeq (b);

        // We declare an iterator on a given sequence.
        Kmer<span>::ModelDirect::Iterator itKmer (model);

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We set the data from which we want to extract kmers.
            itKmer.setData (itSeq->getData());

            // We iterate the kmers.
            for (itKmer.first(); !itKmer.isDone(); itKmer.next())
            {
                cout << model.toString (itKmer->value()) << endl;
                nbKmers++;
            }

            //  We increase the sequences counter.
            nbSequences++;
        }

        // We dump some information about the iterations
        cout << "FOUND " << nbKmers << " kmers in " << nbSequences << " sequences" << endl;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
