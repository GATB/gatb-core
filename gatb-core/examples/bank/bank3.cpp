//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank iteration with a functor                             */
/*                                                                              */
/* This snippet shows how to iterate sequences from a FASTA through a functor.  */
/* Note that this approach is necessary when one wants to parallelize the       */
/* iteration (see snippets on multithreading).                                  */
/*                                                                              */
/********************************************************************************/

// We define a functor that will be called for every iterated sequence.
// The advantages are:
//     - we have a code that focuses on the treatment to do on the sequence;
//       this may be interesting when the corresponding code is long and is
//       likely to be moved in an independent method.
//     - such a functor can be reused in other contexts.

struct Functor {  void operator ()  (Sequence& s)  const
{
    // We dump the data size and the comment
    std::cout << "[" << s.getDataSize() << "] " << s.getComment()  << std::endl;

    // We dump the data
    std::cout << s.toString () << std::endl;
}};

/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a Bank instance defined by a list of filenames
    BankFasta b (argc-1, argv+1);

    // We create an iterator over this bank.
    BankFasta::Iterator it (b);

    // We loop over sequences in a "push" fashion (a functor is called for each sequence)
    it.iterate (Functor());
}
//! [snippet1]
