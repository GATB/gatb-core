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

class DSK
{
public:

    /** */
    struct Info
    {
        Info() : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0) {}

        Info& operator+=  (const Info& other)
        {
            nbSeq         += other.nbSeq;
            dataSeq       += other.dataSeq;
            nbKmers       += other.nbKmers;
            checksumKmers += other.checksumKmers;
            return *this;
        }

        u_int64_t nbSeq;
        u_int64_t dataSeq;
        u_int64_t nbKmers;
        kmer_type checksumKmers;
    };

    /********************************************************************************/
    DSK (const string& filename, size_t kmerSize, size_t nbCores)
        : _bank (filename), _filename(filename),
          _model(kmerSize), _kmerSize(kmerSize),
          _dispatcher(nbCores),
          _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
          _max_disk_space(0), _max_memory(9800), _volume(0), _nb_passes(0), _nb_partitions(0)
    {}

    /********************************************************************************/
    void configure ()
    {
        // We get some estimations about the bank
        _bank.estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

        // We may have to build the binary bank if not already existing.
        buildBankBinary ();

        // We get the available space (in MBytes) of the current directory.
        u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;
        cout << "AVALAIBLE SPACE in current directory '" << System::file().getCurrentDirectory() << "' is " << available_space << " MBytes" << endl;

        _max_disk_space = std::min (available_space/2, _bank.getSize() / 1024 / 1024);
        cout << "MAX DISK SPACE is " << _max_disk_space << " KBytes" << endl;

        if (_max_disk_space == 0)  { _max_disk_space = 10000; }

        // We estimate the number of iterations
        _volume    = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1)) * sizeof(kmer_type) / 1024 / 1024;  // in MBytes
        _nb_passes = ( _volume / _max_disk_space ) + 1;

        cout << "NB PASSES is " << _nb_passes << endl;

        size_t max_open_files = System::file().getMaxFilesNumber() / 2;
        u_int64_t volume_per_pass;

        do
        {
            volume_per_pass = _volume / _nb_passes;
            _nb_partitions  = ( volume_per_pass / _max_memory ) + 1;

            if (_nb_partitions >= max_open_files)   { _nb_passes++;  }
            else                                    { break;         }

        } while (1);
    }

    /********************************************************************************/
    void execute ()
    {
        // We use the provided listener if any
        IteratorListener* progress = new Progress (_estimateSeqNb, "DSK");
        LOCAL (progress);

        // We need an iterator on the FASTA bank.
        Iterator<Sequence>* itBank = _bank.iterator();
        LOCAL (itBank);

        // We declare two kmer iterators for the two banks and a paired one that links them.
        KmerModel::Iterator itKmer (_model);

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<Sequence> itSeq (*itBank, 5*1000);

        // We add a listener to the sequences iterator.
        itSeq.addObserver (*progress);

        u_int64_t total_nbKmers       = 0;
        kmer_type total_checksumKmers = 0;

        // We get some information about the kmers.
        Info infos[_nb_passes];

        /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
        for (size_t p=0; p<_nb_passes; p++)
        {
            /** We update the iterator listener with some information message. */
            progress->setMessage ("DSK %d/%d", p+1, _nb_passes);

            /** We create some functors that will do the actual job. */
            vector<ProcessSequence> functors (_dispatcher.getExecutionUnitsNumber());
            for (size_t i=0; i<functors.size(); i++)   {  functors[i] = ProcessSequence (&_model, _nb_passes, p);  }

            // We get current time stamp
            ITime::Value t0 = System::time().getTimeStamp();

            /** We launch the iteration of the sequences iterator with the created functors. */
            _dispatcher.iterate (itSeq, functors);

            // We get current time stamp
            ITime::Value t1 = System::time().getTimeStamp();

            /** We update the information got during the functors execution. */
            for (size_t i=0; i<functors.size(); i++)   {  infos[p] += functors[i].info;   }

            // We dump some information about the iterated kmers;
            cout << "FOUND " << infos[p].nbKmers << " kmers  for " << infos[p].nbSeq << " sequences  (" << infos[p].dataSeq << " bytes)  "
                << "in " << (t1-t0) << " msec (rate " << (double)infos[p].nbKmers / (double) (t1>t0 ? t1-t0 : 1) << " kmers/msec),  "
                << "checksum is " << infos[p].checksumKmers << " (" << Integer::getName()  << ")"
                << endl;
        }

        /** We get the total number of kmers and the global checksum. */
        for (size_t i=0; i<_nb_passes; i++)
        {
            total_nbKmers       += infos[i].nbKmers;
            total_checksumKmers += infos[i].checksumKmers;
        }

        cout << endl;
        cout << "TOTAL KMERS " << total_nbKmers << "  WITH CHECKSUM " << hex << total_checksumKmers << "  WITH " << Integer::getName()  <<  endl;

        for (size_t p=0; p<_nb_passes; p++)
        {
            cout << "   [" << p << "]  " << 100.0 * (double)infos[p].nbKmers / (double) total_nbKmers << endl;
        }
    }

    /********************************************************************************/
    struct ProcessSequence
    {
        void operator() (Sequence& sequence)
        {
            Hash hash;

            info.nbSeq++;
            info.dataSeq += sequence.getData().size();

            /** We build the kmers from the current sequence. */
            model->build (sequence.getData(), kmers);

            /** We loop over the kmers. */
            for (size_t i=0; i<kmers.size(); i++)
            {
                /** We hash the current kmer. */
                kmer_type h = hash (kmers[i]);

                /** We check whether this kmer has to be processed during the current pass. */
                if ((h % nbPass) != pass)  { continue; }

                info.nbKmers ++;
                info.checksumKmers += h;
            }
        }

        ProcessSequence (KmerModel* model, size_t nbPasses, size_t currentPass)
            : pass(currentPass), nbPass(nbPasses), model(model) {}

        ProcessSequence () : pass(0), nbPass(0), model(0) {}

        vector<kmer_type> kmers;

        Info info;

        size_t pass;
        size_t nbPass;
        KmerModel* model;
    };


private:

    Bank      _bank;
    string    _filename;
    KmerModel _model;
    size_t    _kmerSize;

    ParallelCommandDispatcher _dispatcher;

    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;

    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;

    /** */
    void buildBankBinary ()
    {
        string filenameBin = _filename + ".bin";

        if (System::file().doesExist(filenameBin) == false)
        {
            // We declare a binary bank
            BankBinary bankBin (filenameBin);

            // We declare some job listener.
            Progress progress (_estimateSeqNb, "FASTA to binary conversion");

            // We convert the FASTA bank in binary format
            IProperties* props = BankHelper::singleton().convert (_bank, bankBin, &progress);
            LOCAL (props);
        }
    }
};


/********************************************************************************/

int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"  << endl;
        cerr << "   2) FASTA  bank" << endl;
        cerr << "   3) number of cores" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);

    // We get the number of cores to be used (0 means all cores)
    size_t nbCores = argc >=4 ? atoi(argv[3]) : 0;

    /** We create an instance of DSK class. */
    DSK dsk (filename, kmerSize, nbCores);

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We configure dsk. */
        dsk.configure ();

        /** We execute dsk. */
        dsk.execute ();
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
