#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <iostream>

#include <gatb/tools/math/Integer.hpp>

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

/********************************************************************************/

struct Hash  {   kmer_type operator () (kmer_type& lkmer)
{
    kmer_type kmer_hash;

    kmer_hash  = lkmer ^ (lkmer >> 14);
    kmer_hash  = (~kmer_hash) + (kmer_hash << 18);
    kmer_hash ^= (kmer_hash >> 31);
    kmer_hash  = kmer_hash * 21;
    kmer_hash ^= (kmer_hash >> 11);
    kmer_hash += (kmer_hash << 6);
    kmer_hash ^= (kmer_hash >> 22);

    return kmer_hash;
}};

/********************************************************************************/

struct FunctorIter1
{
    void operator() (Sequence& sequence)
    {
        Hash hash;

        nbSeq++;
        dataSeq += sequence.getData().size();

        /** We build the kmers from the current sequence. */
        model->build (sequence.getData(), kmers);

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            /** We hash the current kmer. */
            kmer_type h = hash (kmers[i]);

            /** We check whether this kmer has to be processed during the current pass. */
            if ((h % nbPass) != pass)  { continue; }

            nbKmers ++;
            checksumKmers = checksumKmers + h;
        }
    }

    FunctorIter1 () : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0), pass(0), nbPass(0), model(0) {}

    vector<kmer_type> kmers;

    u_int64_t nbSeq;
    u_int64_t dataSeq;
    u_int64_t nbKmers;
    kmer_type checksumKmers;

    size_t pass;
    size_t nbPass;
    KmerModel* model;
};

/********************************************************************************/

void iter1 (
    IBank& bank,
    KmerModel& model,
    ICommandDispatcher& dispatcher,
    size_t nbPasses
)
{
    // We use the provided listener if any
    IteratorListener* progress = new Progress (bank.estimateNbSequences(), "Iterating 1");
    LOCAL (progress);

    // We need an iterator on the FASTA bank.
    Iterator<Sequence>* itBank = bank.iterator();
    LOCAL (itBank);

    // We declare two kmer iterators for the two banks and a paired one that links them.
    KmerModel::Iterator itKmer (model);

    // We create an iterator over the paired iterator on sequences
    SubjectIterator<Sequence> itSeq (*itBank, 5*1000);

    if (progress)  {  itSeq.addObserver (*progress);  }

    u_int64_t total_nbKmers       = 0;
    kmer_type total_checksumKmers = 0;

    // We get some information about the kmers.
    u_int64_t nbSeq[nbPasses];
    u_int64_t dataSeq[nbPasses];
    u_int64_t nbKmers[nbPasses];
    kmer_type checksumKmers[nbPasses];

    for (size_t p=0; p<nbPasses; p++)
    {
        /** We update the iterator listener with some information message. */
        progress->setMessage ("DSK %d/%d", p+1, nbPasses);

        // We get current time stamp
        ITime::Value t0 = System::time().getTimeStamp();

        vector<FunctorIter1> functors (dispatcher.getExecutionUnitsNumber());

        for (size_t i=0; i<functors.size(); i++)
        {
            functors[i].pass   = p;
            functors[i].nbPass = nbPasses;
            functors[i].model  = &model;
        }

        nbSeq[p]         = 0;
        dataSeq[p]       = 0;
        nbKmers[p]       = 0;
        checksumKmers[p] = 0;

        /** We launch the iteration. */
        dispatcher.iterate (itSeq, functors);

        for (size_t i=0; i<functors.size(); i++)
        {
            nbSeq[p]         += functors[i].nbSeq;
            dataSeq[p]       += functors[i].dataSeq;
            nbKmers[p]       += functors[i].nbKmers;
            checksumKmers[p] = checksumKmers[p] + functors[i].checksumKmers;

            total_nbKmers       += functors[i].nbKmers;
            total_checksumKmers = total_checksumKmers + functors[i].checksumKmers;
        }

        // We get current time stamp
        ITime::Value t1 = System::time().getTimeStamp();

        // We dump some information about the iterated kmers;
        cout << "FOUND " << nbKmers[p] << " kmers  for " << nbSeq[p] << " sequences  (" << dataSeq[p] << " bytes)  "
            << "in " << (t1-t0) << " msec (rate " << (double)nbKmers[p] / (double) (t1>t0 ? t1-t0 : 1) << " kmers/msec),  "
            << "checksum is " << checksumKmers[p] << " (" << Integer::getName()  << ")"
            << endl;
    }
    cout << endl;
    cout << "TOTAL KMERS " << total_nbKmers << "  WITH CHECKSUM " << hex << total_checksumKmers << "  WITH " << Integer::getName()  <<  endl;

    for (size_t p=0; p<nbPasses; p++)
    {
        cout << "   [" << p << "]  " << 100.0 * (double)nbKmers[p] / (double) total_nbKmers << endl;
    }
}

/********************************************************************************/

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"  << endl;
        cerr << "   2) FASTA  bank" << endl;
        cerr << "   3) binary bank" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);
    string filenameBin = argc <=3 ? (filename + ".bin") : argv[3];

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        ParallelCommandDispatcher dispatcher;

        // We declare a kmer model with a given span size.
        KmerModel model (kmerSize);

        // We declare the FASTA bank
        Bank bank (filename);

        // We declare a binary bank
        BankBinary bankBin (filenameBin);

        // We get some estimation about the bank
        u_int64_t seqNb, seqTotalSize, seqMaxSize;      bank.estimate (seqNb, seqTotalSize, seqMaxSize);

        if (System::file().doesExist(filenameBin) == false)
        {
            // We declare some job listener.
            Progress progress (seqNb, "FASTA to binary conversion");

            // We convert the FASTA bank in binary format
            IProperties* props = BankHelper::singleton().convert (bank, bankBin, &progress);
            LOCAL (props);
        }

        // We get the available space (in MBytes) of the current directory.
        u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;
        cout << "AVALAIBLE SPACE in current directory '" << System::file().getCurrentDirectory() << "' is " << available_space << " MBytes" << endl;

        u_int64_t max_disk_space = std::min (available_space/2, bank.getSize() / 1024 / 1024);
        cout << "MAX DISK SPACE is " << max_disk_space << " KBytes" << endl;

        if (max_disk_space == 0)  { max_disk_space = 10000; }

        // We estimate the number of iterations
        u_int64_t volume    = (seqTotalSize - seqNb * (kmerSize-1)) * sizeof(kmer_type) / 1024 / 1024;  // in MBytes
        u_int32_t nb_passes = ( volume / max_disk_space ) + 1;
        cout << "NB PASSES is " << nb_passes << endl;

cout << "-----> volume=" << volume << "  nb_passes=" << nb_passes << "  max_disk_space=" << max_disk_space << endl;

        u_int64_t volume_per_pass;
        u_int32_t nb_partitions;

        size_t max_open_files = System::file().getMaxFilesNumber() / 2;
        size_t max_memory = 9800;  // in MBytes

        do
        {
            volume_per_pass = volume / nb_passes;
            nb_partitions   = ( volume_per_pass / max_memory ) + 1;

            cout << "volume_per_pass=" << volume_per_pass << "  nb_partitions=" << nb_partitions << endl;

            if (nb_partitions >= max_open_files)    { nb_passes++;  }
            else                                    { break;        }

        } while (1);

        cout << "-----> nb_partitions=" << nb_partitions <<  endl;

        Hash hash;

        // TEST 1
        iter1 (bankBin, model, dispatcher, nb_passes);
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
