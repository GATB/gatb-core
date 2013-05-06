//! [snippet1]
// We include what we need for the test
#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
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
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::math;

//typedef u_int64_t kmer_type;
//typedef ttmath::UInt<64> kmer_type;
//typedef LargeInt<64> kmer_type;
typedef Integer<128> kmer_type;

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

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare a kmer model with a given span size.
        Model<kmer_type> model (kmerSize);

        // We declare the FASTA bank
        Bank bank1 (filename);

        // We declare a binary bank
        BankBinary bank2 (filename + ".bin");

        // We need an iterator on the FASTA bank.
        Bank::Iterator itSeq1 (bank1);

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<Sequence> itSeq1Notif (itSeq1, 1000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        ProgressTimer progressConvert (bank1.estimateNbSequences(), "Iterating sequences during binary conversion");
        itSeq1Notif.addObserver (progressConvert);

        for (itSeq1Notif.first(); !itSeq1Notif.isDone(); itSeq1Notif.next())   {  bank2.insert (*itSeq1);  }   bank2.flush ();

        // We declare two kmer iterators for the two banks and a paired one that links them.
        Model<kmer_type>::Iterator itKmer1 (model, KMER_DIRECT);
        Model<kmer_type>::Iterator itKmer2 (model, KMER_DIRECT);
        PairedIterator<kmer_type,kmer_type> itKmer (itKmer1, itKmer2);

        // We loop the two banks with a paired iterator.
        BankBinary::Iterator itSeq2 (bank2);
        PairedIterator<Sequence,Sequence> itSeqPair (itSeq1, itSeq2);

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<pair<Sequence,Sequence> > itSeq (itSeqPair, 1000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        ProgressTimer progressKmers (bank1.estimateNbSequences(), "Iterating sequences during kmers check");
        itSeq.addObserver (progressKmers);

        // We get some information about the kmers.
        u_int64_t nbKmers       = 0;
        kmer_type checksumKmers = 0;

        // We get current time stamp
        ITime::Value t0 = System::time().getTimeStamp();

        // We loop the sequences of the two banks.
        for (itSeq.first(); !itSeq.isDone();  itSeq.next())
        {
            // We set the data from which we want to extract kmers.
            itKmer1.setData (itSeq1->getData());
            itKmer2.setData (itSeq2->getData());

            // We loop the kmers for the two datas.
            for (itKmer.first(); !itKmer.isDone();  itKmer.next())
            {
                // Here, we should have  (itKmer->first == itKmer->second)
                if (itKmer->first != itKmer->second)  {  cerr << "NO MATCHING KMERS..." << endl; }

                nbKmers       += 1;
                checksumKmers = checksumKmers + itKmer->first;
            }
        }

        // We get current time stamp
        ITime::Value t1 = System::time().getTimeStamp();

        // We dump some information about the iterated kmers;
        cout << "FOUND " << nbKmers << " kmers in " << (t1-t0) << " msec,  checksum is " << checksumKmers  << endl;

        // We remove the binary bank
        System::file().remove (filename + ".bin");
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
