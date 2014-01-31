//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

static const char* STR_URI_BANK_IDS = "-bank-ids";

/********************************************************************************/
/*           Iterate a bank whose sequences are each kmer of a model.           */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankFilter");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,     "bank reference",   true));
    parser.push_back (new OptionOneParam (STR_URI_BANK_IDS,  "file holding indexes of bank",   true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We read the list of indexes. */
        set<size_t> indexes;
        FILE* file = fopen (options->getStr(STR_URI_BANK_IDS).c_str(), "r");
        if (file != 0)
        {
            char buffer[128];
            while (fgets (buffer, sizeof(buffer), file))  {  indexes.insert (atoi(buffer));  }
            fclose (file);
        }

        cout << "found " << indexes.size() << " indexes" << endl;

        /** We open the output bank. */
        string outputBankUri = options->getStr(STR_URI_INPUT) + "_" + System::file().getBaseName (options->getStr(STR_URI_BANK_IDS));
        IBank* outputBank = BankRegistery::singleton().getFactory()->createBank (outputBankUri);
        LOCAL (outputBank);

        /** We loop the input bank. */
        IBank* inputBank = BankRegistery::singleton().getFactory()->createBank (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        inputBank->iterate ([&] (Sequence& seq)
        {
            set<size_t>::iterator lookup = indexes.find (seq.getIndex());
            if (lookup != indexes.end())
            {
                /** The index is ok, we put the sequence into the output bank. */
                outputBank->insert (seq);
            }
        });

        /** We flush the output bank. */
        outputBank->flush();
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
