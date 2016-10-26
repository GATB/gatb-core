//! [snippet1]
// GATB-Core Online Tutorial 
/********************************************************************************/
/*                    Bank iteration with progress information                  */
/*                                                                              */
/* This snippet shows how to iterate sequences from a FASTA with some progress  */
/* information. In this example, we provide our own progress manager.           */
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
// We define a functor that will be called during bank parsing
struct ProgressFunctor : public IteratorListener  {  
  void inc (u_int64_t current){
    // produce a '.' on stdout each time the functor is called 
    // (see the creation of the Sequence Iterator, below)
    std::cout << ".";
  }
};

/********************************************************************************/
// START Application
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

  try
  {
    // We declare an input Bank and use it locally
    IBank* inputBank = Bank::open (argv[1]);
    LOCAL (inputBank);

    // Note also that we have to parameterize the SubjectIterator by the 
    // kind of iterated items (i.e. Sequence) and the processing that has 
    // to be done on each iteration (ProgressFunctor).
    SubjectIterator<Sequence> iter (inputBank->iterator(), 100);

    // We create some listener to be notified every 100 iterations and 
    // attach it to the iterator.
    iter.addObserver (new ProgressFunctor());

    // We loop over sequences.
    for (iter.first(); !iter.isDone(); iter.next())
    {
    // Note that we do nothing inside the sequence iterating loop about 
    // the progression management.
    // In other words, we don't "pollute" the code inside this loop by 
    // presentation concerns and we can therefore focus on the job to 
    // be done on the iterated sequences.
    }
  }
  catch (Exception& e)
  {
    std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
  }
}
//! [snippet1]
