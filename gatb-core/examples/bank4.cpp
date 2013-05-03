//! [snippet1]
// We include what we need for the test
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::tools::dp::impl;

// We a define a functor that will be called during bank parsing
struct ProgressFunctor : public IteratorListener  {  void udpate (u_int64_t current)   {  printf (".");  } };

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails
    try
    {
        // We declare a Bank instance defined by a list of filenames
        Bank b (argc-1, argv+1);

        // We create a sequence iterator for the bank
        Bank::Iterator itSeq (b);

        // Note how we create a Sequence iterator with b.iterator() and give it to the
        // SubjectIterator that will add a notification feature to the Sequence iterator.

        // Note also that we have to parameterize the SubjectIterator by the kind of iterated
        // items (Sequence) and the processing that has to be done on each iteration (ProgressFunctor).
        SubjectIterator<Sequence> iter (itSeq, 10);

        //  We create some listener to be notified every 10 iterations and attach it to the iterator.
        ProgressFunctor fct;    iter.addObserver (fct);

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
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
