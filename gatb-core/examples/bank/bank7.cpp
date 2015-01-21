//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

size_t threshold = 500;

// We a define a functor that will be called during iteration for filtering odd items.
// Here, we keep sequences whose data size is greater than 500.
struct FilterFunctor  {  bool operator ()  (Sequence& seq) const  {  return seq.getDataSize() >= threshold; } };

/********************************************************************************/
/*                                Bank filtering                                */
/*                                                                              */
/* This snippet shows how to iterate a bank and filter some sequences through   */
/* a functor.                                                                   */
/* Note: lambda expressions could be used for the functor (in case the used     */
/* compiler supports it)                                                        */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least the size threshold and the bank file path. Arguments are:" << endl;
        cerr << "   1) bank" << endl;
        cerr << "   2) minimum size threshold (default 500)" << endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);

        // We set the threshold with given parameter
        threshold = atoi (argv[2]);

        // We use another iterator for filtering out some sequences.
        FilterIterator<Sequence,FilterFunctor> itSeq (inputBank->iterator(), FilterFunctor());
        
        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We dump the data size and the comment
            cout << "[" << itSeq->getDataSize() << "] " << itSeq->getComment()  << endl;

            // We dump the data
            cout << itSeq->toString() << endl;
        }
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
