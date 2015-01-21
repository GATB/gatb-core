//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                         Bank conversion with some filtering                  */
/*                                                                              */
/* This snippet shows how to copy a bank to another one with some filtering     */
/* criteria. The criteria is described through a functor.                       */
/* Here, we can keep only every 'modulo' sequence that have no 'N' in data.     */
/*                                                                              */
/********************************************************************************/

// We define a functor for our sequence filtering.
struct FilterFunctor
{
    int modulo;
    FilterFunctor (int modulo) : modulo(modulo) {}
    bool operator ()  (Sequence& seq)
    {
        return (seq.getIndex() % modulo == 0)  &&  (seq.toString().find("N") == string::npos);
    }
};

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
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);

        // We declare an output Bank
        BankFasta outputBank (argv[2]);

        // We get the modulo value (1 by default), defined as static to be retrieved in the functor
        int modulo = argc >= 4 ? atoi (argv[3]) : 1;
        if (modulo<=0)  { modulo=1; }

        // We create a filter on the input bank.
        // Note the first argument of the FilterIterator is an iterator of the input bank and the
        // second argument is the filtering functor.
        FilterIterator<Sequence,FilterFunctor>* itFilter = new FilterIterator<Sequence,FilterFunctor>(
            inputBank->iterator(),
            FilterFunctor (modulo)
        );

        // We create an iterator that will provide progress notification.
        // Note the estimated number of items takes into account the modulo
        ProgressIterator<Sequence> itSeq (itFilter, "Filtering bank", inputBank->estimateNbItems()/modulo );

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We insert the current sequence into the output bank.
            outputBank.insert (*itSeq);
        }

        // We make sure that the output bank is flushed correctly.
        outputBank.flush ();
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
