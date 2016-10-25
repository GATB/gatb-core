//! [snippet1]

// We include GATB-Core
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank iteration with a functor                             */
/*                                                                              */
/* This snippet shows how to iterate sequences from a FASTA through a functor.  */
/*                                                                              */
/* As a reminder, a c++ functor is pretty much just a class which defines the   */
/* operator() function. That lets you create objects which "look like" a        */
/* function. GATB-Core uses that concept to handle various operations on        */
/* objects such as sequences or k-mers.                                         */
/*                                                                              */
/* Note that inGATB, this approach is mandatory when one wants to parallelize   */
/* the iteration (see snippets on multithreading).                              */
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
  // We dump the sequence size and the comment
  std::cout << "[" << s.getDataSize() << "] " << s.getComment()  << std::endl;

  // We dump the sequence
  std::cout << s.toString () << std::endl;
}};

/********************************************************************************/
int main (int argc, char* argv[])
{
  // We check that the user provides at least one option: a Fasta/FastQ file.
  // Online GATB-Tutorial: this argument is automatically filled in with an 
  // appropriate file.
  if (argc < 2)
  {
    std::cerr << "Please, provide a sequence file." << std::endl;
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
