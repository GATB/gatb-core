//! [snippet1]
// We include what we need for the test
#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <iostream>
#include <iomanip>

// We use the required packages
using namespace std;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::misc::impl;

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide an input and output banks names:" << endl;
        cerr << "   1) FASTA file path" << endl;
        cerr << "   2) binary file path" << endl;
        cerr << "   3) dump binary output: 0 for no (default), 1 for yes" << endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare an input Bank
        Bank inputBank (argv[1]);

        // We declare an output Bank
        BankBinary outputBank (argv[2]);

        // We create a sequence iterator for the bank
        Bank::Iterator itInput (inputBank);

        // We create an iterator over the input bank and encapsulate it with progress notification.
        SubjectIterator<Sequence> itSeq (itInput, 100000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        ProgressTimer progress (inputBank.estimateNbSequences(), "Converting input file into binary format");

        itSeq.addObserver (progress);

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We insert the current sequence into the output bank.
            outputBank.insert (*itSeq);
        }

        // We make sure that the output bank is flushed correctly.
        outputBank.flush ();

        if (argc >= 4  && atoi(argv[3])==1)
        {
            // We create an iterator on our binary bank
            BankBinary::Iterator itBinary (outputBank);

            // We iterate the sequences whose data is binary coded
            for (itBinary.first(); !itBinary.isDone(); itBinary.next())
            {
                // We display some (dummy) comment
                cout << itBinary->getComment() << endl;

                for (size_t i=0; i< (itBinary->getDataSize()+3)/4; i++)
                {
                    cout << setfill('0') << setw(2) << hex << (int) (itBinary->getData()[i] & 0xFF) << " ";

                    if ((i+1)%32==0)  { cout << endl; }
                }
                cout << endl;
            }
        }
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
