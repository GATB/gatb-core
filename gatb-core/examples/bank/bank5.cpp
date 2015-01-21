//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank iteration with progress information                  */
/*                                                                              */
/* This snippet shows how to iterate sequences from a bank with some progress   */
/* information. In this example, we use some pre defined progress manager.      */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "you must provide a bank." << std::endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails
    try
    {
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);

        // We create a sequence iterator for the bank with progress information
        ProgressIterator<Sequence> iter (*inputBank, "Iterating sequences");

        // We loop over sequences.
        for (iter.first(); !iter.isDone(); iter.next())
        {
            // Note that we do nothing inside the sequence iterating loop
        }
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
