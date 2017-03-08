//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*           Iterate a bank whose sequences are each kmer of a model.           */
/*                                                                              */
/* This snippet shows how iterate all possible kmers of a given size and dump   */
/* them as output. This program can be viewed as a kmer generator.              */
/*                                                                              */
/* Cmd-line: bank11 -kmer-size <size> -out <fasta file>                         */
/*                                                                              */
/* Sample: bank11 -kmer-size 4 -out /tmp/kmer.fa                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("KmerIteration");
    parser.push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size",          true));
    parser.push_back (new OptionOneParam (STR_URI_OUTPUT,  "output fasta bank",  false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We create a kmers bank with the provided kmer size. */
        BankKmers bank (options->getInt(STR_KMER_SIZE));

        /** We create the output bank. */
        string outputName = options->get(STR_URI_OUTPUT) ?
            options->getStr(STR_URI_OUTPUT) :
            Stringify::format ("bank_k%d.fa", options->getInt(STR_KMER_SIZE));

        BankFasta outputBank (outputName);

        /** We create a Sequence iterator. */
        ProgressIterator<Sequence> it (bank);

        /** We iterate the bank. */
        for (it.first(); !it.isDone(); it.next())  {  outputBank.insert (it.item());  }

        /** We flush the output bank. */
        outputBank.flush();
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
