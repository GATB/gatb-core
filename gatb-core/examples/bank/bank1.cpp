//! [snippet1]
// We include what we need for the test
#include <gatb/bank/impl/Bank.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::bank::impl;

int main (int argc, char* argv[])
{
    // We get the file name from the user arguments
    const char* filename = argc >= 2 ? argv[1] : "";

//! [snippet1_bank]
    // We declare a Bank instance.
    Bank b (filename);

    // We create an iterator over this bank.
    Bank::Iterator it (b);

    // We loop over sequences.
    for (it.first(); !it.isDone(); it.next())
    {
        // In the following, see how we access the current sequence information through
        // the -> operator of the iterator

        // We dump the data size and the comment
        cout << "[" << it->getDataSize() << "] " << it->getComment()  << endl;

        // We dump the data
        cout << it->getDataBuffer() << endl;
    }
//! [snippet1_bank]
}
//! [snippet1]
