//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank with exception management                            */
/*                                                                              */
/* This snippet shows how to open a bank and iterate over sequences.            */
/*                                                                              */
/* Here, we use a more a more generic class to achieve that operation:          */
/* Bank, instead of BankFasta.                                                   */
/*                                                                              */
/* Note: we use here a try/catch block in case the bank opening doesn't work.   */
/*                                                                              */
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

      // We dump the sequence size and the comment
      std::cout << "[" << seq.getDataSize() << "] " << seq.getComment()  << std::endl;

      // We dump the sequence
      std::cout << seq.toString() << std::endl;
    }
  }
  catch (Exception& e)
  {
    std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
  }
}
//! [snippet1]
