#if 0
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
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>

#include <iostream>
#include <map>
#include <math.h>

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

const char* PROP_KMER_SIZE  = "-kmer-size";
const char* PROP_DB         = "-db";
const char* PROP_NB_CORES   = "-nb-cores";
const char* PROP_MAX_MEMORY = "-max-memory";
const char* PROP_NKS        = "-nks";
const char* PROP_PREFIX     = "-prefix";
const char* PROP_VERBOSE    = "-verbose";
const char* PROP_STATS_XML  = "-stats-xml";

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

    /********************************************************************************/
    DSK (IProperties* params)
        : _params(0), _stats(0), _bankBinary(0),
          _nks (3), _prefix("dsk"),
          _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
          _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0)
    {
        setParams (params);

        IProperty* prop;

        if ( (prop = (*params)[PROP_KMER_SIZE])  != 0)  {  _kmerSize   = prop->getInt();     }
        if ( (prop = (*params)[PROP_DB])         != 0)  {  _filename   = prop->getValue();   }
        if ( (prop = (*params)[PROP_MAX_MEMORY]) != 0)  {  _max_memory = prop->getInt();     }
        if ( (prop = (*params)[PROP_NKS])        != 0)  {  _nks        = prop->getInt();     }
        if ( (prop = (*params)[PROP_PREFIX])     != 0)  {  _prefix     = prop->getValue();   }

        /** We create the binary bank holding the reads in binary format. */
        _bankBinary = new BankBinary (_filename + ".bin");

        /** We create our dispatcher. */
        prop = (*params)[PROP_NB_CORES];
        _dispatcher = new ParallelCommandDispatcher (prop ? prop->getInt() : 0);

        /** We create a properties object for gathering statistics. */
        setStats (new Properties());
        _stats->add (0, "dsk");

        /** We add the user parameters to the global stats. */
        _stats->add (1, params);
    }

    /********************************************************************************/
    virtual ~DSK ()
    {
        setParams (0);
        setStats  (0);

        delete _bankBinary;
        delete _dispatcher;
    }

    /********************************************************************************/
    void execute ()
    {
        /** We configure dsk. */
        configure ();

        // We create the sequences iterator.
        Iterator<Sequence>* itSeq = createSequenceIterator (new Progress (_estimateSeqNb, "DSK"));
        LOCAL (itSeq);

        // We create the solid kmers bag
        Bag<kmer_type>* solidKmers = createSolidKmersBag ();
        LOCAL (solidKmers);

        /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
        for (size_t pass=0; pass<_nb_passes; pass++)
        {
            /** 1) We fill the partition files. */
            fillPartitions (pass, itSeq);

            /** 2) We fill the kmers solid file from the partition files. */
            fillSolidKmers (solidKmers);
        }

        /** We flush the solid kmers file. */
        solidKmers->flush();

        /** We add the exec time stats to the global statistics. */
        _stats->add (1, _timeInfo.getProperties("time"));
    }

    /********************************************************************************/
    IProperties& getStats ()  const  { return *_stats; }

private:

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

        /** We gather some statistics. */
        _stats->add (1, "config");
        _stats->add (2, "current directory", System::file().getCurrentDirectory());
        _stats->add (2, "available space",   "%ld", available_space);
        _stats->add (2, "bank size",         "%ld", bankSize);
        _stats->add (2, "sequence number",   "%ld", _estimateSeqNb);
        _stats->add (2, "sequence volume",   "%ld", _estimateSeqTotalSize / MBYTE);
        _stats->add (2, "kmers number",      "%ld", kmersNb);
        _stats->add (2, "kmers volume",      "%ld", _volume);
        _stats->add (2, "max disk space",    "%ld", _max_disk_space);
        _stats->add (2, "nb passes",         "%d",  _nb_passes);
        _stats->add (2, "nb partitions",     "%d",  _nb_partitions);
    }

    /********************************************************************************/
    void fillPartitions (size_t pass, Iterator<Sequence>* itSeq)
    {
        TIME_INFO (_timeInfo, "fill partitions");

        /** We create a kmer model. */
        KmerModel model (_kmerSize);

        /** We create the partition files for the current pass. */
        BagFilePartition<kmer_type> partitions (_nb_partitions, getPartitionUri().c_str());

        /** We create a shared synchronizer for the partitions building. */
        ISynchronizer* synchro = System::thread().newSynchronizer();

        /** We launch the iteration of the sequences iterator with the created functors. */
        _dispatcher->iterate (*itSeq, FillPartitions (model, _nb_passes, pass, partitions, synchro));

        /** We cleanup resources. */
        delete synchro;
    }

    /********************************************************************************/
    void fillSolidKmers (Bag<kmer_type>*  solidKmers)
    {
        TIME_INFO (_timeInfo, "fill solid kmers");

        /** We parse each partition file. */
        for (size_t i=0; i<_nb_partitions; i++)
        {
            char filename[128];  snprintf (filename, sizeof(filename), getPartitionUri().c_str(), i);

            IteratorFile<kmer_type> it (filename);
            vector<kmer_type> kmers;

            it.fill (kmers);

#ifdef OMP
            omptl::sort (kmers.begin (), kmers.end ());
#else
            std::sort (kmers.begin (), kmers.end ());
#endif

            u_int32_t max_couv  = 2147483646;
            u_int32_t abundance = 0;
            kmer_type previous_kmer = kmers.front();

            for (vector<kmer_type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
            {
                kmer_type current_kmer = *itKmers;

                if (current_kmer == previous_kmer)  {   abundance++;  }
                else
                {
                    if (abundance >= _nks  && abundance <= max_couv)
                    {
                        solidKmers->insert (previous_kmer);
                    }
                    abundance = 1;
                }
                previous_kmer = current_kmer;
            }

            if (abundance >= _nks && abundance <= max_couv)
            {
                solidKmers->insert (previous_kmer);
            }
        }
    }

    /********************************************************************************/
    virtual Iterator<Sequence>* createSequenceIterator (IteratorListener* progress)
    {
        // We need an iterator on the FASTA bank.
        Iterator<Sequence>* itBank = _bankBinary->iterator();

        // We create an iterator over the paired iterator on sequences
        SubjectIterator<Sequence>* itSeq = new SubjectIterator<Sequence> (itBank, 5*1000);

        // We add a listener to the sequences iterator.
        itSeq->addObserver (progress);

        // We return the created sequence iterator.
        return itSeq;
    }

    /********************************************************************************/
    virtual Bag<kmer_type>* createSolidKmersBag ()
    {
        /** We delete the solid kmers file. */
        System::file().remove ("solids.bin");

        return new BagCache<kmer_type> (new  BagFile<kmer_type> ("solids.bin"), 5*1000);
    }

    /********************************************************************************/
    void buildBankBinary (Bank& bank)
    {
        TIME_INFO (_timeInfo, "bank conversion");

        // We convert the FASTA bank in binary format
        IProperties* props = BankHelper::singleton().convert (bank, *_bankBinary, new Progress (_estimateSeqNb, "FASTA to binary conversion"));
        LOCAL (props);
    }

    /********************************************************************************/
    string getPartitionUri ()  {  return _prefix + "%d";  }

    /********************************************************************************/
    class FillPartitions
    {
    public:

        void operator() (Sequence& sequence)
        {
            vector<kmer_type> kmers;

            Hash hash;

            /** We build the kmers from the current sequence. */
            model.build (sequence.getData(), kmers);

            /** We loop over the kmers. */
            for (size_t i=0; i<kmers.size(); i++)
            {
                /** We hash the current kmer. */
                kmer_type h = hash (kmers[i]);

                /** We check whether this kmer has to be processed during the current pass. */
                if ((h % nbPass) != pass)  { continue; }

                kmer_type reduced_kmer = h / nbPass;

                /** We compute in which partition this kmer falls into. */
                size_t p = reduced_kmer % _cache.size();

                /** We write the kmer into the bag. */
                _cache[p]->insert (kmers[i]);
            }
        }

        FillPartitions (KmerModel& model, size_t nbPasses, size_t currentPass, BagFilePartition<kmer_type>& partition, ISynchronizer* synchro)
            : pass(currentPass), nbPass(nbPasses), _cache(partition,synchro), model(model)
        {
        }

    private:
        size_t pass;
        size_t nbPass;
        BagCachePartition<kmer_type> _cache;
        KmerModel& model;
    };


    IProperties* _params;
    void setParams (IProperties* params)  { SP_SETATTR(params); }

    IProperties* _stats;
    void setStats (IProperties* stats)  { SP_SETATTR(stats); }

    BankBinary* _bankBinary;

    TimeInfo _timeInfo;

    string   _filename;
    size_t   _kmerSize;
    size_t   _nks;
    string   _prefix;

    ICommandDispatcher* _dispatcher;

    u_int64_t _estimateSeqNb;
    u_int64_t _estimateSeqTotalSize;
    u_int64_t _estimateSeqMaxSize;
    u_int64_t _max_disk_space;
    u_int32_t _max_memory;
    u_int64_t _volume;
    u_int32_t _nb_passes;
    u_int32_t _nb_partitions;
};

/********************************************************************************/
/********************************************************************************/

class Debloom
{
private:

    size_t              _kmerSize;

    IProperties*        _props;
    ICommandDispatcher* _dispatcher;
    TimeInfo            _timeInfo;

    IteratorListener*   _progress;
    void setProgress (IteratorListener* progress)  { SP_SETATTR(progress); }

public:

    /** */
    Debloom (IProperties* props) : _props(props), _progress(0), _stats(0)
    {
        _kmerSize = (*props)[PROP_KMER_SIZE]->getInt();

        size_t nbCores = ((*props)[PROP_NB_CORES] ?  (*props)[PROP_NB_CORES]->getInt() : 0);
        _dispatcher = new ParallelCommandDispatcher (nbCores);

        setStats (new Properties());
        _stats->add (0, "stats");
    }

    /** */
    ~Debloom ()
    {
        delete _dispatcher;
        setProgress (0);
        setStats    (0);
    }

    /** */
    void execute ()
    {
        /** We delete the extension file. */
        System::file().remove ("extension.bin");

        u_int64_t solidFileSize = (System::file().getSize("solids.bin") / sizeof (kmer_type));

        SubjectIterator<kmer_type> itKmers (new IteratorFile<kmer_type> ("solids.bin"), 5*1000);

        setProgress (new Progress (solidFileSize, "iterate solid kmers"));

        itKmers.addObserver (_progress);

        /** We create a bloom with inserted solid kmers. */
        Bloom<kmer_type>* bloom = createBloom (itKmers, solidFileSize);

        BagFile<kmer_type>* extensionFile = new BagFile<kmer_type>("extension.bin");
        LOCAL (extensionFile);

        executeBuildKmerExtension (itKmers, *bloom, extensionFile);

        /** We make sure everything is put into the extension file. */
        extensionFile->flush();

        /** We clean up resources. */
        delete bloom;
    }

    Properties getStats ()
    {
        Properties res;

        res.add (0, "debloom");
        res.add (1, _stats);
        res.add (1, _timeInfo.getProperties("time"));

        return res;
    }

private:

    IProperties* _stats;
    void setStats (IProperties* stats)  { SP_SETATTR(stats); }

    /********************************************************************************/
    Bloom<kmer_type>* createBloom (Iterator<kmer_type>& itKmers, u_int64_t solidFileSize)
    {
        TIME_INFO (_timeInfo, "build kmers bloom");

        double lg2 = log(2);
        float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);

        u_int64_t estimatedBloomSize = solidFileSize * NBITS_PER_KMER;

        /** We instantiate the bloom object. */
        Bloom<kmer_type>* bloom = new BloomSynchronized<kmer_type> (estimatedBloomSize, (int)floorf (0.7*NBITS_PER_KMER));

        _progress->setMessage ("build bloom from solid kmers");

        _dispatcher->iterate (itKmers,  BuildKmerBloom (*bloom));

        /** We gather some statistics. */
        _stats->add (1, "bloom");
        _stats->add (2, "nb solid kmers",   "%d",   solidFileSize);
        _stats->add (2, "nb bits per kmer", "%.3f", NBITS_PER_KMER);
        _stats->add (2, "bloom size",       "%d",   estimatedBloomSize / 8LL);

        return bloom;
    }

    /********************************************************************************/
    void executeBuildKmerExtension (Iterator<kmer_type>& itKmers, Bloom<kmer_type>& bloom, Bag<kmer_type>* extensionBag)
    {
        TIME_INFO (_timeInfo, "build kmers extension");

        KmerModel model (_kmerSize);

        ISynchronizer* synchro = System::thread().newSynchronizer();

        _progress->setMessage ("build extension from solid kmers");

        _dispatcher->iterate (itKmers, BuildKmerExtension (model, bloom, extensionBag, synchro));

        delete synchro;

        /** We gather some statistics. */
        //_stats->add (1, "extension");
    }

    /********************************************************************************/
    struct BuildKmerBloom
    {
        void operator() (const kmer_type& kmer)  {  _bloom.insert(kmer); }
        BuildKmerBloom (Bloom<kmer_type>& bloom)  : _bloom(bloom) {}
        Bloom<kmer_type>& _bloom;
    };

    /********************************************************************************/
    struct BuildKmerExtension
    {
        void operator() (const kmer_type& kmer)
        {
            _itNeighbors.setSource (kmer);

            for (_itNeighbors.first(); !_itNeighbors.isDone(); _itNeighbors.next())
            {
                if (_bloom.contains (*_itNeighbors))
                {
                    _extendBag.insert (*_itNeighbors);
                }
            }
        }

        BuildKmerExtension (KmerModel& model, Bloom<kmer_type>& bloom, Bag<kmer_type>* extendBag, ISynchronizer* synchro)
            : _model(model), _bloom(bloom), _extendBag(extendBag, 5*1000, synchro), _itNeighbors(model)  { }

        ~BuildKmerExtension ()  {  _extendBag.flush();  }

        KmerModel&          _model;
        Bloom<kmer_type>&   _bloom;
        BagCache<kmer_type> _extendBag;
        ModelAbstract<kmer_type>::KmerNeighborIterator _itNeighbors;
    };

};


/********************************************************************************/

int main (int argc, char* argv[])
{
    OptionsParser parser;

    parser.add (new OptionOneParam (PROP_KMER_SIZE,   "size of a kmer",                       true));
    parser.add (new OptionOneParam (PROP_DB,          "URI of the bank",                      true));
    parser.add (new OptionOneParam (PROP_NB_CORES,    "number of cores",                      false));
    parser.add (new OptionOneParam (PROP_MAX_MEMORY,  "max memory",                           false));
    parser.add (new OptionOneParam (PROP_NKS,         "abundance threshold for solid kmers",  false));
    parser.add (new OptionOneParam (PROP_PREFIX,      "prefix URI for temporary files",       false));
    parser.add (new OptionNoParam  (PROP_VERBOSE,     "dump exec info to stdout",             false));
    parser.add (new OptionOneParam (PROP_STATS_XML,   "dump exec info into a XML file",       false));

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We parse the command line arguments. */
        parser.parse (argc, argv);

        /** We get the options as a Properties object. */
        IProperties& props = parser.getProperties();

        /** We read properties from the init file (if any). */
        props.add (0, new Properties (System::info().getHomeDirectory() + string ("/.dskrc")));

        /** We create an instance of DSK class. */
        DSK dsk (&props);

        /** We execute dsk. */
        dsk.execute ();

        /** We create an instance of Debloom class. */
        Debloom debloom (&props);

        /** We execute the debloom. */
        debloom.execute ();

        /** We may have to dump execution information to stdout. */
        if (props[PROP_VERBOSE] != 0)
        {
            RawDumpPropertiesVisitor visit;
            dsk.getStats().accept     (&visit);
            debloom.getStats().accept (&visit);
        }

        /** We may have to dump execution information to stdout. */
        if (props[PROP_STATS_XML] != 0)
        {
            XmlDumpPropertiesVisitor visit (props[PROP_STATS_XML]->getValue());
            dsk.getStats().accept     (&visit);
            debloom.getStats().accept (&visit);
        }
    }

    catch (OptionFailure& e)
    {
        if (parser.saw("-h"))    {   parser.displayHelp   (stdout);   }
        else                     {   parser.displayErrors (stdout);   }
        return EXIT_FAILURE;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
#else

int main (int argc, char* argv[])
{
}
#endif


