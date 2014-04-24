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
/* them as output.                                                              */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("KmerIteration");
    parser.push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We create a kmers bank with the provided kmer size. */
        BankKmers bank (options->getInt(STR_KMER_SIZE));

        /** We create a Sequence iterator. */
        Iterator<Sequence>* it = bank.iterator();
        LOCAL (it);

        /** We iterate the bank. */
        for (it->first(); !it->isDone(); it->next())
        {
            cout << "seq: " << it->item().toString() << endl;
        }
    }

    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
