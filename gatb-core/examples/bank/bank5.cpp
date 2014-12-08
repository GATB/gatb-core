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

        // We create a sequence iterator for the bank with progress information
        ProgressIterator<Sequence> iter (b, "Iterating sequences");

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
