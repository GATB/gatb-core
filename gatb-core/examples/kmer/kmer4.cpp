//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <string.h>

// We use the required packages
using namespace std;

static const size_t span = KMER_SPAN(0);

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how to iterate the kmers of the sequences from a bank.    */
/*                                                                              */
/* Cmd-line: kmer4 -in <fasta/q file> -kmer-size <value> [-verbose]             */
/*                                                                              */
/* Sample: kmer4 -in gatb-core/gatb-core/test/db/reads1.fa -kmer-size 11        */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("KmerTest");
    parser.push_back (new OptionOneParam (STR_URI_INPUT, "bank input",    true));
    parser.push_back (new OptionOneParam (STR_KMER_SIZE, "kmer size",     true));
    parser.push_back (new OptionNoParam  (STR_VERBOSE,   "display kmers", false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We define the max size of a data line in the FASTA output file
        size_t kmerSize = options->getInt(STR_KMER_SIZE);

        bool verbose = options->get(STR_VERBOSE) != 0;

        u_int64_t nbSequences = 0;
        u_int64_t nbKmers     = 0;

        // We open the input bank
        IBank* bank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (bank);

        // We declare a kmer model with a given span size.
        Kmer<span>::ModelDirect model (kmerSize);

        // We declare an iterator on a given sequence.
        Kmer<span>::ModelDirect::Iterator itKmer (model);

        // We create an iterator over this bank.
        ProgressIterator<Sequence> itSeq (*bank);

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We set the data from which we want to extract kmers.
            itKmer.setData (itSeq->getData());

            // We iterate the kmers.
            for (itKmer.first(); !itKmer.isDone(); itKmer.next())
            {
                if (verbose)  {  cout << model.toString (itKmer->value()) << endl;  }
                nbKmers++;
            }

            //  We increase the sequences counter.
            nbSequences++;
        }

        // We dump some information about the iterations
        cout << "FOUND " << nbKmers << " kmers in " << nbSequences << " sequences" << endl;
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
