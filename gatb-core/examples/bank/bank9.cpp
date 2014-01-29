//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <iomanip>

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
        cerr << "   2) output URI" << endl;
        cerr << "   3) modulo (ie keep sequence id % modulo)" << endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare an input Bank
        BankFasta inputBank (argv[1]);

        // We declare an output Bank
        BankFasta outputBank (argv[2]);

        // We get the modulo value (10 by default), defined as static to be retrieved in the functor
        static int modulo = argc >= 4 ? atoi (argv[3]) : 10;

        // We create a sequence iterator for the bank
        BankFasta::Iterator* itInput = new BankFasta::Iterator (inputBank);

        // We define a functor for our sequence filtering.
        struct FilterFunctor  {  bool operator ()  (Sequence& seq)  {
            string data = string (seq.getDataBuffer(), seq.getDataSize());
            return seq.getIndex() % modulo == 0 && strstr (data.c_str(), "N") == 0;
        } } filter;

        // We create a filter on the bank
        FilterIterator<Sequence,FilterFunctor>* itFilter = new FilterIterator<Sequence,FilterFunctor>((Iterator<Sequence>*)itInput, filter);

        // We create an iterator over the input bank and encapsulate it with progress notification.
        SubjectIterator<Sequence> itSeq (itFilter, 100000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        itSeq.addObserver (new ProgressTimer (inputBank.estimateNbSequences(), "Filtering bank"));

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We insert the current sequence into the output bank.
            outputBank.insert (*itSeq);
        }

        // We make sure that the output bank is flushed correctly.
        outputBank.flush ();
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
