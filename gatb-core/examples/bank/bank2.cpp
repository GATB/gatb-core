//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank with exception management                            */
/*                                                                              */
/* This snippet shows how to open a bank and iterate its sequences.             */
/* Note: we use here a try/catch block in case the bank opening doesn't work.   */
/*                                                                              */
/* Cmd-line: bank2 <fasta/q file>                                               */
/*                                                                              */
/* Sample: bank2 gatb-core/gatb-core/test/db/reads1.fa                          */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "you must provide a bank." << std::endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);

        // We create an iterator over this bank.
        Iterator<Sequence>* it = inputBank->iterator();
        LOCAL (it);

        // We loop over sequences.
        for (it->first(); !it->isDone(); it->next())
        {
            // Shortcut
            Sequence& seq = it->item();

            // We dump the data size and the comment
            std::cout << "[" << seq.getDataSize() << "] " << seq.getComment()  << std::endl;

            // We dump the data
            std::cout << seq.toString() << std::endl;
        }
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
