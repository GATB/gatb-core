//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

static const char* STR_URI_SEQ_IDS = "-seq-ids";

// We use a filtering functor that knows which sequence indexes have to be kept.
struct FilterFunctor
{
    set<size_t>& indexes;
    FilterFunctor (set<size_t>& indexes) : indexes(indexes) {}
    bool operator ()  (Sequence& seq) const  {  return indexes.find (seq.getIndex()) != indexes.end();  }
};

/********************************************************************************/
/*           Extract sequences from a bank.                                     */
/*                                                                              */
/* This snippet shows how to extract some sequences for a given list of indexes.*/
/* An "index" i means the i-th sequence inside the input file                   */
/* These indexes are read from an input file, one index per line.               */
/* E.g. if the index files contains:
 *    1                                                                         
      5                                                                         */
/* Then read number 1 and number 5 will be returned */
/*                                                                              */
/* Cmd-line: bank12 -in <fasta/q file> -seq-ids <seq-ids file>                  */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankFilter");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,    "bank reference",               true));
    parser.push_back (new OptionOneParam (STR_URI_SEQ_IDS,  "file holding indexes of bank", true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We read the list of indexes. */
        set<size_t> indexes;
        FILE* file = fopen (options->getStr(STR_URI_SEQ_IDS).c_str(), "r");
        if (file != 0)
        {
            char buffer[128];
            while (fgets (buffer, sizeof(buffer), file))  {  indexes.insert (atoi(buffer));  }
            fclose (file);
        }

        cout << "found " << indexes.size() << " indexes" << endl;

        /** We open the output bank. */
        string outputBankUri = options->getStr(STR_URI_INPUT) + "_" + System::file().getBaseName (options->getStr(STR_URI_SEQ_IDS));
        BankFasta outputBank(outputBankUri);

        /** We loop the input bank. */
        IBank* inputBank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        /** We use another iterator for filtering out some sequences. */
        FilterIterator<Sequence,FilterFunctor> itSeq (inputBank->iterator(), FilterFunctor(indexes));

        /** We loop the sequences. */
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            outputBank.insert (itSeq.item());
        }

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
