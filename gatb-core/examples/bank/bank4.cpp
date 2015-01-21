//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                    Bank iteration with progress information                  */
/*                                                                              */
/* This snippet shows how to iterate sequences from a FASTA with some progress  */
/* information. In this example, we provide our own progress manager.           */
/*                                                                              */
/********************************************************************************/

// We a define a functor that will be called during bank parsing
struct ProgressFunctor : public IteratorListener  {  void inc (u_int64_t current)   {  std::cout << ".";  } };

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "you must provide a bank." << std::endl;
        return EXIT_FAILURE;
    }

    // We define a try/catch block in case some method fails
    try
    {
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);

        // Note also that we have to parameterize the SubjectIterator by the kind of iterated
        // items (Sequence) and the processing that has to be done on each iteration (ProgressFunctor).
        SubjectIterator<Sequence> iter (inputBank->iterator(), 10);

        //  We create some listener to be notified every 10 iterations and attach it to the iterator.
        iter.addObserver (new ProgressFunctor());

        // We loop over sequences.
        for (iter.first(); !iter.isDone(); iter.next())
        {
            // Note that we do nothing inside the sequence iterating loop about the progression management.

            // In other words, we don't "pollute" the code inside this loop by presentation concerns and
            // we can therefore focus on the job to be done on the iterated sequences.
        }
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
