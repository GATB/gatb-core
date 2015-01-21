//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank management                                 */
/*                                                                              */
/* This snippet shows how to open a bank and iterate its sequences.             */
/* Some attributes of the iterated Sequence objects are used.                   */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "you must provide a bank." << std::endl;
        return EXIT_FAILURE;
    }

    // We get the file name from the user arguments
    const char* filename = argc >= 2 ? argv[1] : "";

//! [snippet1_bank]
    // We declare an input Bank and use it locally
    IBank* inputBank = Bank::open (argv[1]);
    LOCAL (inputBank);

    // We create an iterator over this bank.
    Iterator<Sequence>* it = inputBank->iterator();
    LOCAL (it);

    // We loop over sequences.
    for (it->first(); !it->isDone(); it->next())
    {
        // In the following, see how we access the current sequence information through
        // the -> operator of the iterator

        // We dump the data size and the comment
        std::cout << "[" << (*it)->getDataSize() << "] " << (*it)->getComment()  << std::endl;

        // We dump the data
        std::cout << (*it)->toString() << std::endl;
    }
//! [snippet1_bank]
}
//! [snippet1]
