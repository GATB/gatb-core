//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <sstream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                         Bank conversion to binary format                     */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide an input and output banks names:" << endl;
        cerr << "   1) input URI"  << endl;
        cerr << "   2) nb sub banks" << endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare an input Bank
        BankFasta inputBank (argv[1]);

        // We declare an output Bank
        int nbBanks = atoi (argv[2]);
        if (nbBanks <= 0)  { throw Exception ("You should provide a positive number of splits"); }

        // We get the basename of the input bank.
        string bankName = System::file().getBaseName(argv[1]);

        // We create the N output banks.
        vector<IBank*> banks (nbBanks);
        for (size_t i=0; i<banks.size(); i++)
        {
            stringstream ss;  ss << bankName << "_" << i << ".fa";
            banks[i] = BankRegistery::singleton().getFactory()->createBank (ss.str());
        }

        // We create a sequence iterator for the bank
        Iterator<Sequence>* itInput = inputBank.iterator();

        // We create an iterator over the input bank and encapsulate it with progress notification.
        SubjectIterator<Sequence> itSeq1 (itInput, 100000);
        itSeq1.addObserver (new ProgressTimer (inputBank.estimateNbSequences(), "count"));

        // We loop over sequences to get the exact number of sequences.
        size_t nbSequences = 0;
        for (itSeq1.first(); !itSeq1.isDone(); itSeq1.next())  { nbSequences ++; }

        // We create an iterator over the input bank and encapsulate it with progress notification.
        SubjectIterator<Sequence> itSeq2 (itInput, 100000);
        itSeq2.addObserver (new ProgressTimer (inputBank.estimateNbSequences(), "split"));

        // We loop over sequences and put them into the correct split bank.
        for (itSeq2.first(); !itSeq2.isDone(); itSeq2.next())
        {
            size_t id = (*itSeq2).getIndex() / (nbSequences / nbBanks);
            banks[id]->insert (*itSeq2);
        }

        // Some clean up
        for (size_t i=0; i<banks.size(); i++)
        {
            banks[i]->flush();
            delete banks[i];
        }

    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
