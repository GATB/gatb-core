//! [snippet1]
// We include what we need for the test
#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::bank::impl;

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare a Bank instance defined by a list of filenames
        Bank b (argc-1, argv+1);

        // We create an iterator over this bank.
        Bank::Iterator it (b);

        // We loop over sequences.
        for (it.first(); !it.isDone(); it.next())
        {
            // We dump the data size and the comment
            cout << "[" << it->getDataSize() << "] " << it->getComment()  << endl;

            // We dump the data
            cout << it->getDataBuffer() << endl;
        }
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
