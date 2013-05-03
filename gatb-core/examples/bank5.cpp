//! [snippet1]
// We include what we need for the test
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::misc::impl;

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails
    try
    {
        // We declare a Bank instance defined by a list of filenames
        Bank b (argc-1, argv+1);

        // We create a sequence iterator for the bank
        Bank::Iterator itSeq (b);

        //  We create a sequence iterator that notifies listeners every N sequences
        SubjectIterator<Sequence> iter (itSeq, 1000);

        // We create some listener to be notified every N iterations and attach it to the iterator.

        // Note that we get an estimation of the number of sequences in the bank and give it to the
        // functor in order to compute a percentage. Since this is a crude estimation, it is likely
        // that the final percentage will not be exactly 100%
        Progress fct (b.estimateNbSequences(), "Iterating sequences");    iter.addObserver (fct);

        // We loop over sequences.
        for (iter.first(); !iter.isDone(); iter.next())
        {
            // Note that we do nothing inside the sequence iterating loop
        }
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
