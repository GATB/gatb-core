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
/* Cmd-line: bank3 <fasta/q file>                                               */
/*                                                                              */
/* Sample: bank3 gatb-core/gatb-core/test/db/reads1.fa                          */
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
    if (argc < 2)
    {
        std::cerr << "you must provide a bank." << std::endl;
        return EXIT_FAILURE;
    }

    // We declare an input Bank and use it locally
    IBank* inputBank = Bank::open (argv[1]);
    LOCAL (inputBank);

    // We create an iterator over this bank.
    Iterator<Sequence>* it = inputBank->iterator();
    LOCAL (it);

    // We loop over sequences in a "push" fashion (a functor is called for each sequence)
    Functor fct;
    it->iterate (fct);
}
//! [snippet1]
