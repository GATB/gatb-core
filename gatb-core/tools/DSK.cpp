#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <iostream>
#include <map>

#include <gatb/tools/math/Integer.hpp>

#include <omptl/omptl_numeric>
#include <omptl/omptl_algorithm>

#include <omp.h>

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

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

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

struct Partition
{
    Partition (size_t nbPartitions) : partitions(nbPartitions), _nbPartitions(nbPartitions)
    {
        /** We delete the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)
        {
            /** Physical delete. */
            System::file().remove (getFilename(i));
        }

        /** We create the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)
        {
            partitions[i] = new BagFile<kmer_type> (getFilename(i));
        }
    }

    ~Partition ()
    {
        /** We delete the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)
        {
            /** Logical delete. */
         //   delete partitions[i];

            /** Physical delete. */
            System::file().remove (getFilename(i));
        }
    }

    /** */
    void foo ()
    {
        /** We delete the partition files. */
        for (size_t i=0; i<_nbPartitions; i++)
        {
            /** Logical delete. */
            delete partitions[i];
        }
    }

    /** */
    Bag<kmer_type>*& operator[] (size_t idx)  { return partitions[idx]; }

    /** */
    size_t size () const  { return partitions.size(); }

private:

    std::string getFilename (size_t idx)
    {
        char filename[128];  snprintf (filename, sizeof(filename), "foo_partition%d.out", idx);
        return string(filename);
    }

    vector <Bag<kmer_type>*> partitions;
    size_t _nbPartitions;
};

/********************************************************************************/

class DSK
{
public:

    /** */
    struct Info
    {
        Info() : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0) {}

        Info (vector<Info>& infos) : nbSeq(0), dataSeq(0), nbKmers(0), checksumKmers(0)
        {
        	for (size_t i=0; i<infos.size(); i++)  {  *this += infos[i];   }
        }

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
    DSK (const string& filename, size_t kmerSize, ICommandDispatcher& dispatcher, size_t maxMemory)
        : _bankBinary(filename + ".bin"), _filename(filename),
          _model(kmerSize), _kmerSize(kmerSize),
          _dispatcher(dispatcher),
          _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
          _max_disk_space(0), _max_memory(maxMemory), _volume(0), _nb_passes(0), _nb_partitions(0)
    {}

    /********************************************************************************/
    void configure ()
    {
        // We create a Bank instance.
        Bank bank (_filename);

        // We get some estimations about the bank
        bank.estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

        // We may have to build the binary bank if not already existing.
        buildBankBinary (bank);

        // We get the available space (in MBytes) of the current directory.
        u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

        u_int64_t bankSize = bank.getSize() / MBYTE;
        u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));

        _volume = kmersNb * sizeof(kmer_type) / MBYTE;  // in MBytes

        _max_disk_space = std::min (available_space/2, bankSize);

        if (_max_disk_space == 0)  { _max_disk_space = 10000; }

        _nb_passes = ( _volume / _max_disk_space ) + 1;

        size_t max_open_files = System::file().getMaxFilesNumber() / 2;
        u_int64_t volume_per_pass;

        do  {
            volume_per_pass = _volume / _nb_passes;
            _nb_partitions  = ( volume_per_pass / _max_memory ) + 1;

            if (_nb_partitions >= max_open_files)   { _nb_passes++;  }
            else                                    { break;         }

        } while (1);

        cout << "CURRENT DIRECTORY '"    << System::file().getCurrentDirectory() << "'" << endl;
        cout << "AVALAIBLE SPACE     : " << available_space << " MBytes" << endl;
        cout << "BANK SIZE           : " << bankSize << " MBytes" << endl;
        cout << "SEQUENCE NUMBER     : " << _estimateSeqNb << endl;
        cout << "SEQUENCE VOLUME     : " << _estimateSeqTotalSize << " Bytes" << endl;
        cout << "KMER NUMBER         : " << kmersNb << endl;
        cout << "KMERS VOLUME        : " << _volume << " MBytes" << endl;
        cout << "MAX DISK SPACE      : " << _max_disk_space << " MBytes" << endl;
        cout << "MAX MEMORY          : " << _max_memory     << " MBytes" << endl;
        cout << "NB PASSES           : " << _nb_passes << endl;
        cout << "NB PARTITIONS       : " << _nb_partitions << endl;
    }

    /********************************************************************************/
    void execute ()
    {
        // We use the provided listener if any
        IteratorListener* progress = new Progress (_estimateSeqNb, "DSK");
        LOCAL (progress);

        // We need an iterator on the FASTA bank.
        Iterator<Sequence>* itBank = _bankBinary.iterator();
        LOCAL (itBank);

        // We declare two kmer iterators for the two banks and a paired one that links them.
        KmerModel::Iterator itKmer (_model);

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<Sequence> itSeq (*itBank, 5*1000);

        // We add a listener to the sequences iterator.
        itSeq.addObserver (*progress);

        // We get some information about the kmers.
        vector<Info> infos(_nb_passes);

        BagFile<kmer_type> solidKmersFile ("solids.bin");
        BagCache<kmer_type> solidKmers (solidKmersFile, 5*1000);

        /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
        for (size_t p=0; p<_nb_passes; p++)
        {
            /** We update the iterator listener with some information message. */
            progress->setMessage ("DSK %d/%d", p+1, _nb_passes);

            /** We create the partition files. */
            Partition  partitions (_nb_partitions);

            /** 1) We fill the partition files. */
            fillPartitions (p, itSeq, partitions, infos[p]);

            partitions.foo ();

            /** 2) We fill the kmers solid file from the partition files. */
            fillSolidKmers (partitions, solidKmers);
        }

        /** We flush the solid kmers file. */
        solidKmers.flush();

        Info globalInfo (infos);
        cout << endl << "TOTAL KMERS " << globalInfo.nbKmers << "  WITH CHECKSUM " << hex << globalInfo.checksumKmers << "  WITH " << Integer::getName()  <<  endl;
        for (size_t p=0; p<_nb_passes; p++)
        {
            cout << "   [" << p << "]  " << 100.0 * (double)infos[p].nbKmers / (double) globalInfo.nbKmers << endl;
        }

        /** We dump some information about the iterated kmers */
        cout << "nbSolidKmers is " << System::file().getSize(solidKmersFile.getName()) / sizeof(kmer_type)  << endl;
    }

    /********************************************************************************/
    void fillPartitions (
        size_t                  p,
        Iterator<Sequence>&     itSeq,
        Partition&              partition,
        Info&                   info
    )
    {
        /** We create a shared synchronizer for the partitions building. */
        ISynchronizer* synchro = System::thread().newSynchronizer();

        /** We create some functors that will do the actual job. */
        vector<ProcessSequence*> functors (_dispatcher.getExecutionUnitsNumber());
        for (size_t i=0; i<functors.size(); i++)
        {
            functors[i] = new ProcessSequence (&_model, _nb_passes, p, partition, synchro);
        }

        // We get current time stamp
        ITime::Value t0 = System::time().getTimeStamp();

        /** We launch the iteration of the sequences iterator with the created functors. */
        _dispatcher.iterate (itSeq, functors);

        // We get current time stamp
        ITime::Value t1 = System::time().getTimeStamp();

        /** We update the information got during the functors execution. */
        for (size_t i=0; i<functors.size(); i++)
        {
            info += functors[i]->info;

            functors[i]->flush();
        }

        // We dump some information about the iterated kmers;
        cout << "FOUND " << info.nbKmers << " kmers  for " << info.nbSeq << " sequences  (" << info.dataSeq << " bytes)  "
            << "in " << (t1-t0) << " msec (rate " << (double)info.nbKmers / (double) (t1>t0 ? t1-t0 : 1) << " kmers/msec),  "
            << "checksum is " << info.checksumKmers << " (" << Integer::getName()  << ")"
            << endl;

        /** Cleanup. */
        for (size_t i=0; i<functors.size(); i++)  {  delete functors[i]; }

        /** We cleanup resources. */
        if (synchro) { delete synchro; }
    }

    /********************************************************************************/
    void fillSolidKmers (
        Partition&       partitions,
        Bag<kmer_type>&  solidKmers
    )
    {
        /** We parse each partition file. */
        for (size_t i=0; i<_nb_partitions; i++)
        {
            char filename[128];  snprintf (filename, sizeof(filename), "foo_partition%d.out", i);

            IteratorFile<kmer_type> it (filename);

            u_int64_t len = System::file().getSize(filename) / sizeof (kmer_type);
            vector<kmer_type> kmers (len);

ITime::Value tt0 = System::time().getTimeStamp();

            size_t k=0;  for (it.first(); !it.isDone(); it.next())  {  kmers[k++] = *it;  }

ITime::Value tt1 = System::time().getTimeStamp();

            omptl::sort (kmers.begin (), kmers.end ());

ITime::Value tt2 = System::time().getTimeStamp();

            u_int32_t nks = 3;
            u_int32_t max_couv  = 2147483646;
            u_int32_t abundance = 0;
            kmer_type previous_kmer = kmers.front();

            for (vector<kmer_type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
            {
                kmer_type current_kmer = *itKmers;

                if (current_kmer == previous_kmer)  {   abundance++;  }
                else
                {
                    if (abundance >= nks  && abundance <= max_couv)
                    {
                        solidKmers.insert (previous_kmer);
                    }
                    abundance = 1;
                }
                previous_kmer = current_kmer;
            }

            if (abundance >= nks && abundance <= max_couv)
            {
                solidKmers.insert (previous_kmer);
            }
ITime::Value tt3 = System::time().getTimeStamp();

cout << "kmersNb: "  << len << "  " << "fill: " << (tt1-tt0) << "  " << "sort: " << (tt2-tt1) << "  " << "solid: " << (tt3-tt2) << "  " << endl;

        }

    }

    /********************************************************************************/
    class ProcessSequence
    {
    public:

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
                //if ((h & 7 ) != pass)  { continue; }

                kmer_type reduced_kmer = h / nbPass;

                /** We compute in which partition this kmer falls into. */
                size_t p = reduced_kmer % nbPartitions;

                /** We write the kmer into the bag. */
                cache[p]->insert (kmers[i]);

                /** Some statistics. */
                info.nbKmers ++;
                info.checksumKmers += h;
            }
        }

        ProcessSequence (KmerModel* model, size_t nbPasses, size_t currentPass, Partition& partition, ISynchronizer* synchro)
            : pass(currentPass), nbPass(nbPasses), cache(partition.size()), nbPartitions(partition.size()), model(model)
        {
            /** We create the cache for each file bag. */
            for (size_t i=0; i<nbPartitions; i++)
            {
                cache[i] = new BagCache<kmer_type> (*(partition[i]), 5*1000, synchro);
                //cache[i] = new BagCacheSorted<kmer_type> (*(partition[i]), 5*1000, synchro);
            }
        }

        ~ProcessSequence ()
        {
            for (size_t i=0; i<nbPartitions; i++)
            {
                if (cache[i])  {  delete cache[i];  }
            }
            nbPartitions = 0;
        }

        void flush ()
        {
            for (size_t i=0; i<nbPartitions; i++)  {  cache[i]->flush();  }

        }

        ProcessSequence () : pass(0), nbPass(0), cache(1), nbPartitions(1), model(0)
        {
            for (size_t i=0; i<nbPartitions; i++)  {  cache[i] = 0; }
        }

        vector<kmer_type> kmers;

        Info info;

        size_t pass;
        size_t nbPass;
        vector < Bag<kmer_type>* >  cache;
        size_t nbPartitions;
        KmerModel* model;
    };


private:

    BankBinary _bankBinary;

    string    _filename;
    KmerModel _model;
    size_t    _kmerSize;

    ICommandDispatcher& _dispatcher;

    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;

    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;

    /** */
    void buildBankBinary (Bank& bank)
    {
        string filenameBin = _filename + ".bin";

        if (System::file().doesExist(filenameBin) == false)
        {
            // We declare a binary bank
            BankBinary bankBin (filenameBin);

            // We declare some job listener.
            Progress progress (_estimateSeqNb, "FASTA to binary conversion");

            // We convert the FASTA bank in binary format
            IProperties* props = BankHelper::singleton().convert (bank, bankBin, &progress);
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
        cerr << "   4) memory (in MBytes)" << endl;
        return EXIT_FAILURE;
    }

    // We define the max size of a data line in the FASTA output file
    size_t kmerSize = atoi(argv[1]);

    // We get the URI of the FASTA bank
    string filename (argv[2]);

    // We get the number of cores to be used (0 means all cores)
    size_t nbCores = argc >=4 ? atoi(argv[3]) : 0;

    // We get the memory to use (in MBytes)
    size_t maxMemory = argc >=5 ? atoi(argv[4]) : 0;

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We create a command dispatcher. */
        ICommandDispatcher* dispatcher = 0;
        if (nbCores==1)  { dispatcher =  new SerialCommandDispatcher   ();         }
        else             { dispatcher =  new ParallelCommandDispatcher (nbCores);  }
        LOCAL (dispatcher);

    	/** We create an instance of DSK class. */
    	DSK dsk (filename, kmerSize, *dispatcher, maxMemory);

        /** We configure dsk. */
        dsk.configure ();

        /** We execute dsk. */
        dsk.execute ();
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    cout << "finished..." << endl;

    return EXIT_SUCCESS;
}
