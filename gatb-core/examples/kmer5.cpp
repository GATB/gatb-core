//! [snippet1]
// We include what we need for the test
#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <iostream>
#include <string.h>

// We use the required packages
using namespace std;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

// We need a functor that links two iterators; it updates the inner loop iterator (on kmers)
// with the data of the outer iterator (on sequences).
struct Update { void operator() (Iterator<kmer_type>* itKmer, Sequence* seq)
{
    // We have to recover the real type of the iterator (lost due to genericity of InnerIterator)
    static_cast<KmerModel::Iterator*> (itKmer)->setData (seq->getData());
}};

// We need a functor that shows some iteration progress
struct ProgressFunctor : public IteratorListener
{
    ProgressFunctor (u_int64_t n) : total(n>0?n:1) {}
    ~ProgressFunctor () { cout << endl; }
    void inc  (u_int64_t current)   {  cout << current << " " << 100.0*(float)current / (float)total << "\r";  cout.flush(); }
    u_int64_t total;
};

////////////////////////////////////////////////////////////////////////////////
int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        u_int64_t nbKmers = 0;

        // We declare a Bank instance defined by a list of filenames
        Bank bank (argc-1, argv+1);

        // We declare a kmer model with a given span size.
        KmerModel model (27);

        // We create a sequence iterator for the bank
        Bank::Iterator* itBank = new Bank::Iterator (bank);

        // We create a sequences iterator and makes it launch progression notif every N iterations.
        SubjectIterator<Sequence> itSeq (itBank, 100);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        itSeq.addObserver (new ProgressFunctor (bank.estimateNbSequences()));

        // We declare a kmer iterator for the model
        KmerModel::Iterator itKmer (model);

        // We create a compound iterator that iterates kmer from sequences
        CompoundIterator<Sequence,kmer_type,Update> it (itSeq, itKmer, Update());

        // We iterate the kmers.
        for (it.first(); !it.isDone(); it.next())   {  nbKmers++;  }

        // We dump some information about the iterations
        cout << "FOUND " << nbKmers << " kmers" << endl;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
