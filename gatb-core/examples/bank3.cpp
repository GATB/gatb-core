//! [snippet1]
// We include what we need for the test
#include <gatb/bank/impl/Bank.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

// We define a functor that will be called for every iterated sequence.
// The advantages are:
//     - we have a code that focuses on the treatment to do on the sequence;
//       this may be interesting when the corresponding code is long and is
//       likely to be moved in an independant method.
//     - such a functor can be reused in other contexts.

struct Functor {  void operator ()  (Sequence& s)  const
{
    // We dump the data size and the comment
    cout << "[" << s.getDataSize() << "] " << s.getComment()  << endl;

    // We dump the data
    cout << s.getDataBuffer () << endl;
}};

int main (int argc, char* argv[])
{
    // We declare a Bank instance defined by a list of filenames
    Bank b (argc-1, argv+1);

    // We create an iterator over this bank.
    Bank::Iterator it (b);

    // We loop over sequences in a "push" fashion (a functor is called for each sequence)
    Functor fct;    it.iterate (fct);
}
//! [snippet1]
