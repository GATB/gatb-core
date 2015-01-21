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
        cerr << "   1) bank file path" << endl;
        cerr << "   2) binary file path" << endl;
        cerr << "   3) dump binary output: 0 for no (default), 1 for yes" << endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        //! [snippet8_binary]
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);

        // We declare an output Bank
        BankBinary outputBank (argv[2]);

        // We create a sequence iterator on the input bank (with progress information).
        ProgressIterator<Sequence> itSeq (*inputBank, "Converting input file into binary format");

        // We insert each sequence of the input bank into the output bank.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())  {  outputBank.insert (itSeq.item());  }

        // We make sure that the output bank is flushed correctly.
        outputBank.flush ();
        //! [snippet8_binary]

        if (argc >= 4  && atoi(argv[3])==1)
        {
            // We create an iterator on our binary bank
            Iterator<Sequence>* itBinary = outputBank.iterator();
            LOCAL (itBinary);

            // We iterate the sequences whose data is binary coded
            for (itBinary->first(); !itBinary->isDone(); itBinary->next())
            {
                // Shortcut
                Sequence& seq = itBinary->item();

                // We display some (dummy) comment
                cout << seq.getComment() << endl;

                // In binary format, each byte holds 4 nucleotides.

                for (size_t i=0; i< (seq.getDataSize()+3)/4; i++)
                {
                    cout << "[" << setfill('0') << setw(2) << hex << (int) (seq.getData()[i] & 0xFF) << " ";

                    // We convert back from binary to ASCII
                    char b = (seq.getData()[i] & 0xFF) ;

                    static char tableASCII[] = {'A', 'C', 'T', 'G'};

                    cout << tableASCII[(b>>6)&3]
                         << tableASCII[(b>>4)&3]
                         << tableASCII[(b>>2)&3]
                         << tableASCII[(b>>0)&3]
                         << "] ";

                    if ((i+1)%16==0)  { cout << endl; }
                }
                cout << endl;
            }
        }
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
