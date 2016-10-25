//! [snippet1]

// We include GATB-Core
#include <gatb/gatb_core.hpp>
#include <iostream>

// We use the required packages
using namespace std;

// sequence size threshold; see its use in the FilterFunctor
size_t threshold = 500;

// We a define a functor that will be called during iteration for 
// filtering sequences.
struct FilterFunctor  {  
  bool operator ()  (Sequence& seq) const  {  
    // Here, we keep sequences whose size is greater than 500.
    return seq.getDataSize() >= threshold; 
  }
};

/********************************************************************************/
/*                                Bank filtering                                */
/*                                                                              */
/* This snippet shows how to iterate a bank and filter some sequences through   */
/* a functor.                                                                   */
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

  try
  {
    // We declare an input Bank and use it locally
    IBank* inputBank = Bank::open (argv[1]);
    LOCAL (inputBank);

    // We use another iterator for filtering out some sequences.
    FilterIterator<Sequence,FilterFunctor> itSeq (inputBank->iterator(), FilterFunctor());

    // We loop over sequences.
    for (itSeq.first(); !itSeq.isDone(); itSeq.next())
    {
      // We dump the sequence size and its comment
      cout << "[" << itSeq->getDataSize() << "] " << itSeq->getComment()  << endl;

      // We could also dump the sequence
      //cout << itSeq->toString() << endl;
    }
  }
  catch (Exception& e)
  {
    cerr << "EXCEPTION: " << e.getMessage() << endl;
  }
}
//! [snippet1]
