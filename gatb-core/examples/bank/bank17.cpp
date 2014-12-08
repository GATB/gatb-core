//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Bank album                                      */
/*                                                                              */
/* A BankAlbum file is a bank defined by a list of URI of other banks.          */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We get the file name from the user arguments
    const char* filename = argc >= 2 ? argv[1] : "";

    //! [snippet17_album]
    // We declare a Bank instance.
    BankAlbum bank (filename);

    // We dump some information about the bank
    std::cout << "cummulated files sizes : " << bank.getSize() << std::endl;

    // We create an iterator on the bank
    Iterator<Sequence>* it = bank.iterator();
    LOCAL (it);

    // We iterate the sequences of the bank
    for (it->first(); !it->isDone(); it->next())
    {
        // We dump some information about the sequence.
        std::cout << "comment " << it->item().getComment() << std::endl;
    }
    //! [snippet17_album]
}
//! [snippet1]
