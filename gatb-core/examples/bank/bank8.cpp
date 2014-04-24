//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <iomanip>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                         Bank conversion to binary format                     */
/*                                                                              */
/* This snippet shows how to convert a bank from a format to another one. The   */
/* main idea is to iterate the input bank and to insert each sequence into the  */
/* output bank. We use progress information to get feedback about the iteration */
/* progression.                                                                 */
/*                                                                              */
/********************************************************************************/
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
        BankFasta inputBank (argv[1]);

        // We declare an output Bank
        BankBinary outputBank (argv[2]);

        // We create a sequence iterator for the bank
        BankFasta::Iterator* itInput = new BankFasta::Iterator (inputBank);

        // We create an iterator over the input bank and encapsulate it with progress notification.
        SubjectIterator<Sequence> itSeq (itInput, 100000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        itSeq.addObserver (new ProgressTimer (inputBank.estimateNbItems(), "Converting input file into binary format"));

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
