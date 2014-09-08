/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

//#define PROTO_COMP
#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <cmath>

#define DEBUG(a)  //printf a

/********************************************************************************/
// We use the required packages
/********************************************************************************/
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::kmer::impl;


/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

/********************************************************************************/
static const char* progressFormat1 = "DSK: Pass %d/%d, Step 1: partitioning    ";
static const char* progressFormat2 = "DSK: Pass %d/%d, Step 2: counting kmers  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :  
R: why are we repeating those long initializations for all constructors?
R: I think it has to do with: http://stackoverflow.com/questions/761917/handling-a-class-with-a-long-initialization-list-and-multiple-constructors
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>::SortingCountAlgorithm ()
    : Algorithm("dsk", 0, 0),
      _storage(0),
      _bank(0),
      _kmerSize(0), _abundance(0),
      _partitionType(0), _nbCores(0), _prefix(""),
      _progress (0),
      _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
      _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
      _histogram (0), _histogramUri(""),
      _partitionsStorage(0), _partitions(0), _totalKmerNb(0), _solidKmers(0)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>::SortingCountAlgorithm (
    Storage* storage,
    gatb::core::bank::IBank* bank,
    size_t      kmerSize,
    size_t      abundance,
    u_int32_t   max_memory,
    u_int64_t   max_disk_space,
    size_t      nbCores,
    size_t      partitionType,
    const std::string& prefix,
    const std::string& histogramUri,
    gatb::core::tools::misc::IProperties* options
)
  : Algorithm("dsk", nbCores, options),
    _storage(storage),
    _bank(0),
    _kmerSize(kmerSize), _abundance(abundance),
    _partitionType(partitionType), _nbCores(nbCores), _prefix(prefix),
    _progress (0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(max_disk_space), _max_memory(max_memory), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
    _histogram (0), _histogramUri(histogramUri),
    _partitionsStorage(0), _partitions(0), _totalKmerNb(0)
{
    setBank (bank);

    /** We create the collection corresponding to the solid kmers output. */
    setSolidKmers (& (*_storage)("dsk").getCollection<Count> ("solid"));

    /** We set the histogram instance. */
    setHistogram (new Histogram  (10000, & (*_storage)("dsk").getCollection<Histogram::Entry>("histogram") ));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>::SortingCountAlgorithm (tools::storage::impl::Storage& storage)
  : Algorithm("dsk", 0, 0),
    _storage(0),
    _bank(0),
    _kmerSize(0), _abundance(0),
    _partitionType(0), _nbCores(0), _prefix(""),
    _progress (0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
    _histogram (0), _histogramUri(""),
    _partitionsStorage(0), _partitions(0), _totalKmerNb(0), _solidKmers(0)
{
    Group& group = storage(this->getName());

    /** We create the collection corresponding to the solid kmers output. */
    setSolidKmers (& group.getCollection<Count> ("solid"));

    string xmlString = group.getProperty ("xml");
    stringstream ss; ss << xmlString;   getInfo()->readXML (ss);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>::~SortingCountAlgorithm ()
{
    setProgress          (0);
    setBank              (0);
    setPartitionsStorage (0);
    setPartitions        (0);
    setSolidKmers        (0);
    setHistogram         (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>& SortingCountAlgorithm<span>::operator= (const SortingCountAlgorithm& s)
{
    if (this != &s)
    {
        _storage                = s._storage;
        _kmerSize               = s._kmerSize;
        _abundance              = s._abundance;
        _partitionType          = s._partitionType;
        _nbCores                = s._nbCores;
        _prefix                 = s._prefix;
        _estimateSeqNb          = s._estimateSeqNb;
        _estimateSeqTotalSize   = s._estimateSeqTotalSize;
        _estimateSeqMaxSize     = s._estimateSeqMaxSize;
        _max_disk_space         = s._max_disk_space;
        _max_memory             = s._max_memory;
        _volume                 = s._volume;
        _nb_passes              = s._nb_passes;
        _nb_partitions          = s._nb_partitions;
        _current_pass           = s._current_pass;
        _histogramUri           = s._histogramUri;
        _totalKmerNb            = s._totalKmerNb;

        setBank                 (s._bank);
        setProgress             (s._progress);
        setHistogram            (s._histogram);
        setPartitionsStorage    (s._partitionsStorage);
        setPartitions           (s._partitions);
        setSolidKmers           (s._solidKmers);
    }
    return *this;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void SortingCountAlgorithm<span>::execute ()
{
    /** We retrieve the actual number of cores. */
    _nbCores = getDispatcher()->getExecutionUnitsNumber();
    assert (_nbCores > 0);

    /** We configure dsk by computing the number of passes and partitions we will have
     * according to the allowed disk and memory space. */
    configure (_bank);

    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = _bank->iterator();
    LOCAL (itSeq);

    /** We configure the progress bar. */
    setProgress ( createIteratorListener (2 * _volume * MBYTE / sizeof(Type), "counting kmers"));
    _progress->init ();

    /*************************************************************/
    /*                         MAIN LOOP                         */
    /*************************************************************/
    /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
    for (_current_pass=0; _current_pass<_nb_passes; _current_pass++)
    {
        DEBUG (("SortingCountAlgorithm<span>::execute  pass [%ld,%d] \n", _current_pass+1, _nb_passes));

        /** 1) We fill the partition files. */
        fillPartitions (_current_pass, itSeq);

        /** 2) We fill the kmers solid file from the partition files. */
        fillSolidKmers (_solidKmers->bag());
    }

    _progress->finish ();

    /** We flush the solid kmers file. */
    _solidKmers->bag()->flush();

    /** We save the histogram if any. */
    _histogram->save ();
	
	
    /** compute auto cutoff **/
    _histogram->compute_threshold ();

	/** store auto cutoff and corresponding number of solid kmers **/
	Collection<NativeInt64>& storecutoff =   (*_storage)("dsk").getCollection<NativeInt64>("cutoff") ;
	storecutoff.insert(_histogram->get_solid_cutoff());
	storecutoff.flush();
	
	Collection<NativeInt64>& storesolids =   (*_storage)("dsk").getCollection<NativeInt64>("nbsolidsforcutoff") ;
	storesolids.insert(_histogram->get_nbsolids_auto());
	storesolids.flush();
	


    /** We want to remove physically the partitions. */
    _partitions->remove ();

    u_int64_t nbSolids = _solidKmers->iterable()->getNbItems();

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "kmers_nb_valid",     "%ld", _totalKmerNb);
    getInfo()->add (2, "kmers_nb_solid",     "%ld", nbSolids);
    getInfo()->add (2, "kmers_nb_weak",      "%ld", _totalKmerNb - nbSolids);
    if (_totalKmerNb > 0)  {  getInfo()->add (2, "kmers_percent_weak", "%.1f", 100.0 - 100.0 * (double)nbSolids / (double)_totalKmerNb  );  }

    getInfo()->add (1, getTimeInfo().getProperties("time"));

    /** We save (as metadata) some information. */
    (*_storage)("dsk").addProperty ("kmer_size", Stringify::format("%d", _kmerSize));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void SortingCountAlgorithm<span>::configure (IBank* bank)
{
    float load_factor = 0.7;

    // optimism == 1 mean that we guarantee worst case the memory usage,
    // any value above assumes that, on average, a k-mer will be seen 'optimism' times
    int optimism = 0; // guarantees always work, above 0 : assumes kmer seen (optimism+1) times (with optim 1  'OAHash: max rehashes..'  happened sometimes)

    /** We get some information about the bank. */
    bank->estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));
    u_int64_t bankSize = _estimateSeqTotalSize / MBYTE;

    _volume = kmersNb * sizeof(Type) / MBYTE;  // in MBytes
    if (_volume == 0)   { _volume = 1; }    // tiny files fix

    if (_max_disk_space == 0)  { _max_disk_space = std::min (available_space/2, bankSize);  }
    if (_max_disk_space == 0)  { _max_disk_space = 10000; }

    if (_max_memory == 0)  {  _max_memory = System::info().getMemoryProject(); }
    if (_max_memory == 0)  {  _max_memory = 1000; }

    assert (_max_disk_space > 0);
    _nb_passes = ( _volume / _max_disk_space ) + 1;

    size_t max_open_files = System::file().getMaxFilesNumber() / 2;
    u_int64_t volume_per_pass;

    do  {
        assert (_nb_passes > 0);
        volume_per_pass = _volume / _nb_passes;

        assert (_max_memory > 0);
        _nb_partitions  = ( (volume_per_pass*_nbCores) / _max_memory ) + 1;

        if (_partitionType == 0)
        {
            _nb_partitions = (u_int32_t) ceil((float) _nb_partitions / load_factor);
            _nb_partitions = ((_nb_partitions * OAHash<Type>::size_entry()) + sizeof(Type)-1) / sizeof(Type); // also adjust for hash overhead
            _nb_partitions = std::max ((_nb_partitions/(optimism+1)), (u_int32_t)1);
        }

        if (_nb_partitions >= max_open_files)   { _nb_passes++;  }
        else                                    { break;         }

    } while (1);

    /** We gather some statistics. */
    getInfo()->add (1, "config");
    getInfo()->add (2, "kmer_size",         "%ld", _kmerSize);
    getInfo()->add (2, "abundance",         "%ld", _abundance);
    getInfo()->add (2, "available_space",   "%ld", available_space);
    getInfo()->add (2, "bank_size",         "%ld", bankSize);
    getInfo()->add (2, "sequence_number",   "%ld", _estimateSeqNb);
    getInfo()->add (2, "sequence_volume",   "%ld", _estimateSeqTotalSize / MBYTE);
    getInfo()->add (2, "kmers_number",      "%ld", kmersNb);
    getInfo()->add (2, "kmers_volume",      "%ld", _volume);
    getInfo()->add (2, "max_disk_space",    "%ld", _max_disk_space);
    getInfo()->add (2, "max_memory",        "%ld", _max_memory);
    getInfo()->add (2, "nb_passes",         "%d",  _nb_passes);
    getInfo()->add (2, "nb_partitions",     "%d",  _nb_partitions);
    getInfo()->add (2, "nb_bits_per_kmer",  "%d",  Type::getSize());
    getInfo()->add (2, "nb_cores",          "%d",  getDispatcher()->getExecutionUnitsNumber());
    getInfo()->add (2, "partition_type",    "%d",  _partitionType);
}

/********************************************************************************/

template<size_t span>
class FillPartitions
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type                  Type;
    typedef typename Kmer<span>::ModelCanonical        Model;
    typedef typename Kmer<span>::ModelCanonical::Kmer  Kmer;

    void operator() (Sequence& sequence)
    {
        /** We build the kmers from the current sequence. */
        if (model.build (sequence.getData(), kmers) == false)  { return; }

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            /** We hash the current kmer. */
            u_int64_t h = oahash (kmers[i].value());

            /** We check whether this kmer has to be processed during the current pass. */
            if ((h % nbPass) != pass)  { continue; }

            u_int64_t reduced_kmer = h / nbPass;

            /** We compute in which partition this kmer falls into. */
            size_t p = reduced_kmer % nbPartitions;

            /** We write the kmer into the bag. */
            _partition[p].insert (kmers[i].value());

            nbWrittenKmers++;
        }

        if (nbWrittenKmers > 500000)   {  _progress.inc (nbWrittenKmers);  nbWrittenKmers = 0;  }
    }

    FillPartitions (Model& model, size_t nbPasses, size_t currentPass, Partition<Type>* partition, u_int32_t max_memory, IteratorListener* progress)
        : model(model), pass(currentPass), nbPass(nbPasses), nbPartitions(partition->size()), nbWrittenKmers(0),
#ifdef PROTO_COMP
      _partition (*partition,1<<12,max_memory, 0),
#else
      _partition (*partition,1<<12,0),
#endif
          _progress  (progress,System::thread().newSynchronizer())  {}

private:

    /** Local resources. */
    Model&    model;
    size_t    pass;
    size_t    nbPass;
    size_t    nbPartitions;
    size_t    nbWrittenKmers;
    Data      binaryData;
    vector<Kmer> kmers;

    /** Shared resources (must support concurrent accesses). */ //PartitionCacheSorted
#ifdef PROTO_COMP
    PartitionCacheSorted<Type> _partition;
#else
    PartitionCache<Type> _partition;
#endif
    
    ProgressSynchro _progress;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void SortingCountAlgorithm<span>::fillPartitions (size_t pass, Iterator<Sequence>* itSeq)
{
    TIME_INFO (getTimeInfo(), "fill_partitions");

    DEBUG (("SortingCountAlgorithm<span>::fillPartitions  pass \n", pass));

    /** We create a kmer model. */
    Model model (_kmerSize);

    /** We delete the previous partitions storage. */
    if (_partitionsStorage)  { _partitionsStorage->remove (); }

    /** We create the partition files for the current pass. */
#ifdef PROTO_COMP
    setPartitionsStorage (StorageFactory(STORAGE_COMPRESSED_FILE).create ("partitions", true, false));
    //setPartitionsStorage (StorageFactory(STORAGE_GZFILE).create ("partitions", true, false));
#else
    setPartitionsStorage (StorageFactory(STORAGE_FILE).create ("partitions", true, false));
#endif
    
	setPartitions        (0); // close the partitions first, otherwise new files are opened before  closing parti from previous pass

    setPartitions        ( & (*_partitionsStorage)().getPartition<Type> ("parts", _nb_partitions));

    /** We update the message of the progress bar. */
    _progress->setMessage (progressFormat1, _current_pass+1, _nb_passes);

    /** We launch the iteration of the sequences iterator with the created functors. */
    // R to E: typical question from a puzzled developer: what's this last argument of iterate() (15*1000)? so to answer this, I need to find the definition of iterate(), where is it? maybe to answer this, where does getDispatcher come from? any way we could search http://gatb-core.gforge.inria.fr/ to quickly find out?
    getDispatcher()->iterate (itSeq, FillPartitions<span> (model, _nb_passes, pass, _partitions, _max_memory, _progress), 15*1000);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
std::vector<size_t> SortingCountAlgorithm<span>::getNbCoresList ()
{
    std::vector<size_t> result;

    for (size_t p=0; p<_nb_partitions; )
    {
        size_t i=0;  for (i=0; i<_nbCores && p<_nb_partitions; i++, p++)  {}
        result.push_back (i);
    }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void SortingCountAlgorithm<span>::fillSolidKmers (Bag<Count>*  solidKmers)
{
    TIME_INFO (getTimeInfo(), "fill_solid_kmers");

    DEBUG (("SortingCountAlgorithm<span>::fillSolidKmers\n"));

    /** We update the message of the progress bar. */
    _progress->setMessage (progressFormat2, _current_pass+1, _nb_passes);

    /** We retrieve the list of cores number for dispatching N partitions in N threads.
     *  We need to know these numbers for allocating the N maps according to the maximum allowed memory.
     */
    vector<size_t> coreList = getNbCoresList();

    size_t p = 0;
    for (size_t i=0; i<coreList.size(); i++)
    {
        vector<ICommand*> cmds;

        size_t currentNbCores = coreList[i];
        assert (currentNbCores > 0);

        /** We correct the number of memory per map according to the max allowed memory.
         * Note that _max_memory has initially been divided by the user provided cores number. */
        u_int64_t mem = (_max_memory*MBYTE)/currentNbCores;

        ISynchronizer* synchro = System::thread().newSynchronizer();
        LOCAL (synchro);

        for (size_t j=0; j<currentNbCores; j++, p++)
        {
            ICommand* cmd = 0;

            /** We get the length of the current partition file. */
            size_t partitionLen = (*_partitions)[p].getNbItems(); // hmm this will be unknown for a gz file, maybe estimation possible du coup utilise le hash

            /* Get the memory taken by this partition if loaded for sorting */
            uint64_t memoryPartition = partitionLen * sizeof(Type);

            if (memoryPartition >= mem)
            {
                cmd = new PartitionsByHashCommand<span>   (solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb, _abundance, _progress, mem);
            }
            else
            {
                cmd = new PartitionsByVectorCommand<span> (solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb, _abundance, _progress);
            }

            cmds.push_back (cmd);
        }

        getDispatcher()->dispatchCommands (cmds, 0);
    }
}

/********************************************************************************/

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
