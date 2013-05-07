#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BankKmerIterator.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/math/ttmath/ttmath.h>
#include <iostream>
#include <string.h>

// We use the required packages
using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc::impl;

typedef u_int64_t kmer_type;
//typedef ttmath::UInt<2> kmer_type;
//typedef LargeInt<2> kmer_type;
//typedef Integer<2> kmer_type;
//typedef __uint128_t kmer_type;


/********************************************************************************/

struct HashNull  {   kmer_type operator () (kmer_type& lkmer)
{
    return 0;
}};

struct Hash  {   kmer_type operator () (kmer_type& lkmer)
{
    kmer_type kmer_hash;
    kmer_hash = lkmer ^ (lkmer >> 14);
    kmer_hash = (~kmer_hash) + (kmer_hash << 18);
    kmer_hash = kmer_hash ^ (kmer_hash >> 31);
    kmer_hash = kmer_hash * 21;
    kmer_hash = kmer_hash ^ (kmer_hash >> 11);
    kmer_hash = kmer_hash + (kmer_hash << 6);
    kmer_hash = kmer_hash ^ (kmer_hash >> 22);
    return kmer_hash;
}};

/********************************************************************************/

template <typename Functor> void iter1 (IBank& bank, Model<kmer_type>& model, Functor& hash, IteratorListener* progress=0)
{
    // We use the provided listener if any
    LOCAL (progress);

    // We need an iterator on the FASTA bank.
    Iterator<Sequence>* itBank = bank.iterator();
    LOCAL (itBank);

    // We declare two kmer iterators for the two banks and a paired one that links them.
    Model<kmer_type>::Iterator itKmer (model, KMER_MINIMUM);

    // We get some information about the kmers.
    u_int64_t nbKmers       = 0;
    kmer_type checksumKmers = 0;

#if 1
    // We create an iterator over the paired iterator on sequences
    SubjectIterator<Sequence> itSeq (*itBank, 10*1000);

    if (progress)  {  itSeq.addObserver (*progress);  }
#else
    Iterator<Sequence>& itSeq = *itBank;
#endif

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    for (itSeq.first(); !itSeq.isDone(); itSeq.next())
    {
        // We set the data from which we want to extract kmers.
        itKmer.setData (itSeq->getData());

        // We loop the kmers for the two datas.
        for (itKmer.first(); !itKmer.isDone();  itKmer.next())
        {
            nbKmers       += 1;
            checksumKmers = checksumKmers + hash (*itKmer);
        }
    }

    // We get current time stamp
    ITime::Value t1 = System::time().getTimeStamp();

    // We dump some information about the iterated kmers;
    cout << "FOUND " << nbKmers << " kmers "
         << "in " << (t1-t0) << " msec (rate " << (double)nbKmers / (double) (t1-t0) << " kmers/msec),  "
         << "checksum is " << checksumKmers
         << endl;
}

/********************************************************************************/

template <typename Functor> void iter2 (IBank& bank, Model<kmer_type>& model, Functor& hash, IteratorListener* progress=0)
{
    // We use the provided listener if any
    LOCAL (progress);

    BankKmerIterator<kmer_type> itKmerBank (bank, model, KMER_MINIMUM);

    // We get some information about the kmers.
    u_int64_t nbKmers       = 0;
    kmer_type checksumKmers = 0;

    // We get current time stamp
    ITime::Value t0 = System::time().getTimeStamp();

    if (progress)  {  itKmerBank.addObserver (*progress);  }

    for (itKmerBank.first(); !itKmerBank.isDone();  itKmerBank.next())
    {
        nbKmers      += 1;
        checksumKmers = checksumKmers + hash (*itKmerBank);
    }

    // We get current time stamp
    ITime::Value t1 = System::time().getTimeStamp();

    // We dump some information about the iterated kmers;
    cout << "FOUND " << nbKmers << " kmers "
         << "in " << (t1-t0) << " msec (rate " << (double)nbKmers / (double) (t1-t0) << " kmers/msec),  "
         << "checksum is " << checksumKmers
         << endl;
}

/********************************************************************************/

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"  << endl;
        cerr << "   2) FASTA bank" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);
    string filenameBin = filename + ".bin";

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare a kmer model with a given span size.
        Model<kmer_type> model (kmerSize);

        // We declare the FASTA bank
        Bank bank (filename);

        // We declare a binary bank
        BankBinary bankBin (filenameBin);

        if (System::file().doesExist(filenameBin) == false)
        {
            // We declare some job listener.
            Progress progress (bank.estimateNbSequences(), "FASTA to binary conversion");

            // We convert the FASTA bank in binary format
            IProperties* props = BankHelper::singleton().convert (bank, bankBin, &progress);
            LOCAL (props);
        }

        //HashNull hash;
        Hash hash;

        // TEST 1
        iter1 (bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 1"));

        // TEST 2
        iter2 (bankBin, model, hash, new Progress (bank.estimateNbSequences(), "Iterating 2"));

        // We remove the binary bank
        //System::file().remove (filename + ".bin");
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
