//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <string.h>

// We use the required packages
using namespace std;

typedef LargeInt<1> kmer_type;

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
        Iterator<Sequence>* itSeq1 = bank1.iterator();

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<Sequence> itSeq1Notif (itSeq1, 1000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        itSeq1Notif.addObserver (new ProgressTimer (bank1.estimateNbSequences(), "Iterating sequences during binary conversion"));

        for (itSeq1Notif.first(); !itSeq1Notif.isDone(); itSeq1Notif.next())   {  bank2.insert (itSeq1->item());  }   bank2.flush ();

        // We declare two kmer iterators for the two banks and a paired one that links them.
        Model<kmer_type>::Iterator itKmer1 (model);
        Model<kmer_type>::Iterator itKmer2 (model);
        PairedIterator<kmer_type,kmer_type> itKmer (itKmer1, itKmer2);

        // We loop the two banks with a paired iterator.
        Iterator<Sequence>* itSeq2 = bank2.iterator();
        LOCAL (itSeq2);

        PairedIterator<Sequence,Sequence>* itSeqPair = new PairedIterator<Sequence,Sequence> (*itSeq1, *itSeq2);
        LOCAL (itSeqPair);

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<pair<Sequence,Sequence> > itSeq (itSeqPair, 1000);

        // We create some listener to be notified every N iterations and attach it to the iterator.
        itSeq.addObserver (new ProgressTimer (bank1.estimateNbSequences(), "Iterating sequences during kmers check"));

        // We get some information about the kmers.
        u_int64_t nbKmers       = 0;
        kmer_type checksumKmers = 0;

        // We get current time stamp
        ITime::Value t0 = System::time().getTimeStamp();

        // We loop the sequences of the two banks.
        for (itSeq.first(); !itSeq.isDone();  itSeq.next())
        {
            // We set the data from which we want to extract kmers.
            itKmer1.setData ((*itSeq1)->getData());
            itKmer2.setData ((*itSeq2)->getData());

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
        cout << "FOUND " << nbKmers << " kmers in " << (t1-t0) << " msec,  checksum is " << checksumKmers  <<  "  ("
             << kmer_type::getName() << ")" << endl;

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
