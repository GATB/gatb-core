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
    // We define a try/catch block in case some method fails
    try
    {
        // We declare a Bank instance defined by a list of filenames
        BankFasta b (argc-1, argv+1);

        // We create a sequence iterator for the bank
        BankFasta::Iterator* itSeq = new BankFasta::Iterator (b);

        // Note how we create a Sequence iterator with b.iterator() and give it to the
        // SubjectIterator that will add a notification feature to the Sequence iterator.

        // Note also that we have to parameterize the SubjectIterator by the kind of iterated
        // items (Sequence) and the processing that has to be done on each iteration (ProgressFunctor).
        SubjectIterator<Sequence> iter (itSeq, 10);

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
    catch (gatb::core::system::Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
}
//! [snippet1]
