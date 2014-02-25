/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/math/Integer.hpp>

#include <gatb/bank/impl/Banks.hpp>

#include <math.h>
#include <algorithm>

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
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>::SortingCountAlgorithm ()
    : Algorithm("dsk", 0, 0),
      _storage(0),
      _bank(0),
      _kmerSize(0), _nks(0),
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
    size_t      nks,
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
    _kmerSize(kmerSize), _nks(nks),
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
        _nks                    = s._nks;
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
    int optimism = 1;

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
    getInfo()->add (2, "nks",               "%ld", _nks);
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
    typedef typename Kmer<span>::Model Model;
    typedef typename Kmer<span>::Type  Type;

    void operator() (Sequence& sequence)
    {
        /** We build the kmers from the current sequence. */
        if (model.build (sequence.getData(), kmers) == false)  { return; }

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            /** We hash the current kmer. */
            Type h = oahash (kmers[i]);

            /** We check whether this kmer has to be processed during the current pass. */
            if ((h % nbPass) != pass)  { continue; }

            Type reduced_kmer = h / nbPass;

            /** We compute in which partition this kmer falls into. */
            size_t p = reduced_kmer % nbPartitions;

            /** We write the kmer into the bag. */
            _partition[p].insert (kmers[i]);

            nbWrittenKmers++;
        }

        if (nbWrittenKmers > 500000)   {  _progress.inc (nbWrittenKmers);  nbWrittenKmers = 0;  }
    }

    FillPartitions (Model& model, size_t nbPasses, size_t currentPass, Partition<Type>* partition, IteratorListener* progress)
        : model(model), pass(currentPass), nbPass(nbPasses), nbPartitions(partition->size()), nbWrittenKmers(0),
          _partition (*partition,1<<12,System::thread().newSynchronizer()),
          _progress  (progress,System::thread().newSynchronizer())  {}

private:

    /** Local resources. */
    Model&    model;
    size_t    pass;
    size_t    nbPass;
    size_t    nbPartitions;
    size_t    nbWrittenKmers;
    Data      binaryData;
    vector<Type> kmers;

    /** Shared resources (must support concurrent accesses). */
    PartitionCache<Type> _partition;
    ProgressSynchro      _progress;
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
    setPartitionsStorage (StorageFactory(STORAGE_FILE).createStorage ("partitions", true, false));
    setPartitions        ( & (*_partitionsStorage)().getPartition<Type> ("parts", _nb_partitions));

    /** We update the message of the progress bar. */
    _progress->setMessage (progressFormat1, _current_pass+1, _nb_passes);

    /** We launch the iteration of the sequences iterator with the created functors. */
    getDispatcher()->iterate (itSeq, FillPartitions<span> (model, _nb_passes, pass, _partitions, _progress), 15*1000);
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
class PartitionsCommand : public ICommand, public system::SmartPointer
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    PartitionsCommand (
        SortingCountAlgorithm<span>& algo,
        Bag<Count>* solidKmers,
        Iterable<Type>&   partition,
        IHistogram*    histogram,
        ISynchronizer* synchro,
        u_int64_t&     totalKmerNbRef
    )
        : _nks(algo.getNks()),
          _solidKmers(solidKmers, 10*1000, synchro),
          _partition(partition),
          _histogram(histogram),
          _progress(algo.getProgress(), synchro),
          _totalKmerNb(0),
          _totalKmerNbRef(totalKmerNbRef) {}

    ~PartitionsCommand()  { _totalKmerNbRef += _totalKmerNb; }

protected:
    size_t              _nks;
    BagCache<Count>      _solidKmers;
    Iterable<Type>&      _partition;
    HistogramCache      _histogram;
    ProgressSynchro     _progress;
    u_int64_t           _totalKmerNb;
    u_int64_t&          _totalKmerNbRef;

    void insert (const Count& kmer)
    {
        u_int32_t max_couv  = 2147483646;

        _totalKmerNb++;

        /** We should update the abundance histogram*/
        _histogram.inc (kmer.abundance);

        /** We check that the current abundance is in the correct range. */
        if (kmer.abundance >= this->_nks  && kmer.abundance <= max_couv)  {  this->_solidKmers.insert (kmer);  }
    }
};

/********************************************************************************/
/** */
template<size_t span>
class PartitionsByHashCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    PartitionsByHashCommand (
        SortingCountAlgorithm<span>& algo,
        Bag<Count>*      solidKmers,
        Iterable<Type>&  partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        u_int64_t       hashMemory
    )
        : PartitionsCommand<span> (algo, solidKmers, partition, histogram, synchro, totalKmerNbRef), _hashMemory(hashMemory)  {}

    void execute ()
    {
        size_t count=0;

        /** We need a map for storing part of solid kmers. */
        OAHash<Type> hash (_hashMemory);

        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL(it);

        for (it->first(); !it->isDone(); it->next())
        {
            hash.increment (it->item());

            /** Some display. */
            if (++count == 100000)  {  this->_progress.inc (count);  count=0; }
        }

        /** We loop over the solid kmers map. */
        Iterator < Abundance<Type> >* itKmerAbundance = hash.iterator();
        LOCAL (itKmerAbundance);

        for (itKmerAbundance->first(); !itKmerAbundance->isDone(); itKmerAbundance->next())
        {
            /** We may add this kmer to the solid kmers bag. */
           this->insert ((Count&) itKmerAbundance->item());
        }
    }

private:
    u_int64_t _hashMemory;
};

/********************************************************************************/
/** */
template<size_t span>
class PartitionsByVectorCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    PartitionsByVectorCommand (
        SortingCountAlgorithm<span>& algo,
        Bag<Count>*  solidKmers,
        Iterable<Type>&    partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef
    )
        : PartitionsCommand<span> (algo, solidKmers, partition, histogram, synchro, totalKmerNbRef)
          {}

    void execute ()
    {
        /** We get the length of the current partition file. */
        size_t partitionLen = this->_partition.getNbItems();

        /** We check that we got something. */
        if (partitionLen == 0)  { throw Exception ("DSK: no solid kmers found"); }

        /** We resize our vector that will be filled with the partition file content.
         * NOTE: we add an extra item and we will set it to the maximum kmer value. */
        kmers.resize (1 + partitionLen);

        /** We directly fill the vector from the current partition file. */
        Iterator<Type>* it = this->_partition.iterator();  LOCAL (it);
        size_t idx = 0;
        for (it->first(); !it->isDone(); it->next(), idx++) { kmers[idx] = it->item(); }

        /** We set the extra item to a max value, so we are sure it will sorted at the last location.
         * This trick allows to avoid extra treatment after the loop that computes the kmers abundance. */
        kmers[partitionLen] = ~0;

        /** We sort the vector. */
        std::sort (kmers.begin (), kmers.end ());

        u_int32_t abundance = 0;
        Type previous_kmer = kmers.front();

        /** We loop over the sorted solid kmers. */
        for (typename vector<Type>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
        {
            if (*itKmers == previous_kmer)  {   abundance++;  }
            else
            {
                this->insert (Count (previous_kmer, abundance) );

                abundance     = 1;
                previous_kmer = *itKmers;
            }
        }

        /** We update the progress bar. */
        this->_progress.inc (kmers.size());
    }

private:

    vector<Type> kmers;
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

            if (_partitionType == 0)
            {
                cmd = new PartitionsByHashCommand<span>   (*this, solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb, mem);
            }
            else
            {
                cmd = new PartitionsByVectorCommand<span> (*this, solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb);
            }

            cmds.push_back (cmd);
        }

        getDispatcher()->dispatchCommands (cmds, 0);
    }
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class SortingCountAlgorithm <32>;
template class SortingCountAlgorithm <64>;
template class SortingCountAlgorithm <96>;
template class SortingCountAlgorithm <128>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
