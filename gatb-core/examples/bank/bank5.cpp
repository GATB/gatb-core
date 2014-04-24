//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank iteration with progress information                  */
/*                                                                              */
/* This snippet shows how to iterate sequences from a FASTA with some progress  */
/* information. In this example, we use some pre defined progress manager.      */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails
    try
    {
        // We declare a Bank instance defined by a list of filenames
        BankFasta b (argc-1, argv+1);

        // We create a sequence iterator for the bank
        BankFasta::Iterator* itSeq = new BankFasta::Iterator (b);

        //  We create a sequence iterator that notifies listeners every N sequences
        SubjectIterator<Sequence> iter (itSeq, 1000);

        // We create some listener to be notified every N iterations and attach it to the iterator.

        // Note that we get an estimation of the number of sequences in the bank and give it to the
        // functor in order to compute a percentage. Since this is a crude estimation, it is likely
        // that the final percentage will not be exactly 100%
        iter.addObserver (new ProgressTimer (b.estimateNbItems(), "Iterating sequences"));

        // We loop over sequences.
        for (iter.first(); !iter.isDone(); iter.next())
        {
            // Note that we do nothing inside the sequence iterating loop
        }
    }
    catch (gatb::core::system::Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
