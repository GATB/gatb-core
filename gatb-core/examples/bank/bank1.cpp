//! [snippet1]
// GATB-Core Online Tutorial 
/********************************************************************************/
/*                    Simple Bank management                                    */
/*                                                                              */
/* In GATB-Core, a bank denotes a set of DNA sequences. And these sequences     */
/* are supposed to be in a sequence file: so, we use a Bank object to read      */
/* a sequence file.                                                             */
/*                                                                              */
/* This snippet shows how to open a FASTA file and iterate over sequences.      */
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
// START Application
int main (int argc, char* argv[])
{
  // We get the file name from the user arguments
  // Online GATB-Tutorial: this argument is automatically filled in with an 
  // appropriate file.
  const char* filename = argc >= 2 ? argv[1] : "";

  //! [snippet1_bank]
  // We declare a BankFasta instance.
  BankFasta b (filename);

  // We create an iterator over this bank. So, a bank denotes a set of
  // sequences contained in a file.
  BankFasta::Iterator it (b);

  // We loop over sequences.
  for (it.first(); !it.isDone(); it.next())
  {
    // In the following, see how we access the current sequence information through
    // the -> operator of the iterator

    // We dump the data size and the comment
    std::cout << "[" << it->getDataSize() << "] " << it->getComment()  << std::endl;

    // We dump the data
    std::cout << it->toString() << std::endl;
  }
  //! [snippet1_bank]
}
//! [snippet1]
