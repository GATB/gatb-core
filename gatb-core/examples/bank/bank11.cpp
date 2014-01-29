//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*           Iterate a bank whose sequences are each kmer of a model.           */
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

        /** We iterate the bank. */
        bank.iterate ([&] (Sequence& seq)
        {
            cout << "seq: " << seq.toString() << endl;
        });
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
