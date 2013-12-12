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
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "you must provide at least the size threshold and the FASTA file path. Arguments are:" << endl;
        cerr << "   1) minimum size threshold (default 500)" << endl;
        cerr << "   2) FASTA files" << endl;
        return EXIT_FAILURE;
    }

    // We set the threshold with given parameter
    threshold = atoi (argv[1]);

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare a Bank instance defined by a list of filenames
        Bank b (argc-2, argv+2);

        // We declare a functor for filtering items.
        FilterFunctor filter;

        // We use another iterator for filtering out some sequences.
        FilterIterator<Sequence,FilterFunctor> itSeq (b.iterator(), filter);
        
        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We dump the data size and the comment
            cout << "[" << itSeq->getDataSize() << "] " << itSeq->getComment()  << endl;

            // We dump the data
            cout << itSeq->getDataBuffer() << endl;
        }
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
