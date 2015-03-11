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
static const char* progressFormat0 = "DSK: estimating nb distinct kmers        ";
static const char* progressFormat1 = "DSK: Pass %d/%d, Step 1: partitioning    ";
static const char* progressFormat2 = "DSK: Pass %d/%d, Step 2: counting kmers  ";
static const char* progressFormat3 = "DSK: Collecting stats on read sample   ";

static u_int64_t DEFAULT_MINIMIZER = 1000000000 ;

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
      _kmerSize(0), _abundance(make_pair(0,~0)),
      _partitionType(0), _minimizerType(0), _repartitionType(0), _nbCores(0), _prefix(""),
      _progress (0),
      _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
      _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
      _histogram (0),
      _partitionsStorage(0), _partitions(0), _totalKmerNb(0), _solidCounts(0), _solidKmers(0) ,_nbCores_per_partition(1) ,_nb_partitions_in_parallel(0),
      _flagEstimateNbDistinctKmers(false), _estimatedDistinctKmerNb(0), _solidityKind(KMER_SOLIDITY_DEFAULT),
      _min_auto_threshold(3)
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
    std::pair<size_t,size_t> abundance,
    u_int32_t   max_memory,
    u_int64_t   max_disk_space,
    size_t      nbCores,
    KmerSolidityKind solidityKind,
    size_t      histogramMax,
    size_t      partitionType,
    size_t      minimizerType,
    size_t      repartitionType,
    size_t      minimizerSize,
    const std::string& prefix,
    gatb::core::tools::misc::IProperties* options
)
  : Algorithm("dsk", nbCores, options),
    _storage(storage),
    _bank(0),
    _kmerSize(kmerSize), _minim_size(minimizerSize), _abundance(abundance),
    _partitionType(partitionType), _minimizerType(minimizerType), 
    _repartitionType(repartitionType),
    _nbCores(nbCores), _prefix(prefix),
    _progress (0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(max_disk_space), _max_memory(max_memory), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
    _histogram (0),
    _partitionsStorage(0), _partitions(0), _totalKmerNb(0), _solidCounts(0), _solidKmers(0) ,_nbCores_per_partition (1) ,_nb_partitions_in_parallel (nbCores),
    _flagEstimateNbDistinctKmers(false),  _estimatedDistinctKmerNb(0),
    _solidityKind(solidityKind),
    _min_auto_threshold(3)
{
    setBank (bank);

    /** We create the collection corresponding to the solid kmers output. */
	// setSolidCounts (& (*_storage)("dsk").getCollection<Count> ("solid"));

    /** We set the histogram instance. */
    setHistogram (new Histogram  (histogramMax, & (*_storage)("dsk").getCollection<Histogram::Entry>("histogram") ));
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
    _storage(&storage),
    _bank(0),
    _kmerSize(0), _minim_size(0), _abundance(make_pair(0,~0)),
    _partitionType(0), _minimizerType(0), 
    _repartitionType(0),
    _nbCores(0), _prefix(""),
    _progress (0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(0), _max_memory(0), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
    _histogram (0),
    _partitionsStorage(0), _partitions(0), _totalKmerNb(0), _solidCounts(0), _solidKmers(0) ,_nbCores_per_partition(1),_nb_partitions_in_parallel(0),
    _flagEstimateNbDistinctKmers(false),_estimatedDistinctKmerNb(0),
    _solidityKind(KMER_SOLIDITY_DEFAULT),
    _min_auto_threshold(3)
{
    Group& group = (*_storage)(this->getName());

    /** We create the collection corresponding to the solid kmers output. */
    setSolidCounts (& group.getPartition<Count> ("solid"));

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
    setSolidCounts       (0);
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
        _minimizerType          = s._minimizerType;
        _repartitionType        = s._repartitionType;
        _nbCores                = s._nbCores;
        _nbCores_per_partition  = s._nbCores_per_partition;
        _nb_partitions_in_parallel = s._nb_partitions_in_parallel;
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
        _totalKmerNb            = s._totalKmerNb;
        _estimatedDistinctKmerNb   = s._estimatedDistinctKmerNb;
        _flagEstimateNbDistinctKmers = s._flagEstimateNbDistinctKmers;
        _solidityKind           = s._solidityKind;

        setBank                 (s._bank);
        setProgress             (s._progress);
        setHistogram            (s._histogram);
        setPartitionsStorage    (s._partitionsStorage);
        setPartitions           (s._partitions);
        setSolidCounts          (s._solidCounts);
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
IOptionsParser* SortingCountAlgorithm<span>::getOptionsParser (bool mandatory)
{
    IOptionsParser* parser = new OptionsParser ("kmer count");

    parser->push_back (new OptionOneParam (STR_URI_INPUT,         "reads file", mandatory ));
    parser->push_back (new OptionOneParam (STR_KMER_SIZE,         "size of a kmer",                           false,  "31"    ));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN,"min abundance threshold for solid kmers",  false,  "3"     ));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX,"max abundance threshold for solid kmers",  false,  "4294967295"));
    parser->push_back (new OptionOneParam (STR_HISTOGRAM_MAX,     "max number of values in kmers histogram",  false, "10000"));
    parser->push_back (new OptionOneParam (STR_SOLIDITY_KIND,     "way to compute solids (sum, min or max)",  false, "sum"));
    parser->push_back (new OptionOneParam (STR_MAX_MEMORY,        "max memory (in MBytes)",                   false, "2000"));
    parser->push_back (new OptionOneParam (STR_MAX_DISK,          "max disk   (in MBytes)",                   false, "0"));
    parser->push_back (new OptionOneParam (STR_URI_SOLID_KMERS,   "output file for solid kmers",              false));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT,        "output file",                              false));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT_DIR,    "output directory",                         false,  "."));
    parser->push_back (new OptionOneParam (STR_MINIMIZER_TYPE,    "minimizer type (0=lexi, 1=freq)",          false,  "0"));
    parser->push_back (new OptionOneParam (STR_MINIMIZER_SIZE,    "size of a minimizer",                      false,  "8"));
    parser->push_back (new OptionOneParam (STR_REPARTITION_TYPE,  "minimizer repartition (0=unordered, 1=ordered)",false,  "0"));

    return parser;
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
    _nb_partitions_in_parallel = _nbCores;
    assert (_nbCores > 0);

    /** We configure dsk by computing the number of passes and partitions we will have
     * according to the allowed disk and memory space. */
    configure (_bank);

    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = _bank->iterator();
    LOCAL (itSeq);

    /** We configure the progress bar. Note that we create a ProgressSynchro since this progress bar
     * may me modified by several threads at the same time. */
    setProgress (new ProgressSynchro (
        createIteratorListener (2 * _volume * MBYTE / sizeof(Type), "counting kmers"),
        System::thread().newSynchronizer())
    );
    _progress->init ();

    /** We create the PartiInfo instance. */
    PartiInfo<5> pInfo (_nb_partitions, _minim_size);

    /*************************************************************/
    /*                         MAIN LOOP                         */
    /*************************************************************/
    /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
    for (_current_pass=0; _current_pass<_nb_passes; _current_pass++)
    {
        DEBUG (("SortingCountAlgorithm<span>::execute  pass [%ld,%d] \n", _current_pass+1, _nb_passes));

        pInfo.clear();

        /** 1) We fill the partition files. */
        fillPartitions (_current_pass, itSeq, pInfo);

        /** 2) We fill the kmers solid file from the partition files. */
        fillSolidKmers (pInfo);
    }

    _progress->finish ();

    /** We flush the solid kmers file. */
    _solidCounts->flush();

    /** We save the histogram if any. */
    _histogram->save ();

    /** compute auto cutoff **/
    _histogram->compute_threshold (_min_auto_threshold);

    /** store auto cutoff and corresponding number of solid kmers **/
    Collection<NativeInt64>& storecutoff =   (*_storage)("dsk").getCollection<NativeInt64>("cutoff") ;
    storecutoff.insert(_histogram->get_solid_cutoff());
    storecutoff.flush();

    Collection<NativeInt64>& storesolids =   (*_storage)("dsk").getCollection<NativeInt64>("nbsolidsforcutoff") ;
    storesolids.insert(_histogram->get_nbsolids_auto());
    storesolids.flush();

    /** We want to remove physically the partitions. */
    _partitions->remove ();

    u_int64_t nbSolids = _solidCounts->getNbItems();

    /** We gather some statistics. */
    if (_bankStats.sequencesNb > 0)
    {
        getInfo()->add (1, "bank");
        getInfo()->add (2, "bank_uri",          "%s",   _bank->getId().c_str());
        getInfo()->add (2, "bank_size",         "%lld", _bank->getSize());
        getInfo()->add (2, "bank_total_nt",     "%lld", _bankStats.sequencesTotalLength);
        getInfo()->add (2, "sequences");
        getInfo()->add (3, "seq_number",        "%ld",  _bankStats.sequencesNb);
        getInfo()->add (3, "seq_size_min",      "%ld",  _bankStats.sequencesMinLength);
        getInfo()->add (3, "seq_size_max",      "%ld",  _bankStats.sequencesMaxLength);
        getInfo()->add (3, "seq_size_mean",     "%.1f", _bankStats.getSeqMean());
        getInfo()->add (3, "seq_size_deviation","%.1f", _bankStats.getSeqDeviation());
        getInfo()->add (2, "kmers");
        getInfo()->add (3, "kmers_nb_valid",   "%lld", _bankStats.kmersNbValid);
        getInfo()->add (3, "kmers_nb_invalid", "%lld", _bankStats.kmersNbInvalid);
    }

    getInfo()->add (1, "stats");

    getInfo()->add (2, "kmers");
    getInfo()->add (3, "kmers_nb_distinct",  "%ld", _totalKmerNb);
    getInfo()->add (3, "kmers_nb_solid",     "%ld", nbSolids);
    getInfo()->add (3, "kmers_nb_weak",      "%ld", _totalKmerNb - nbSolids);
    if (_totalKmerNb > 0)  {  getInfo()->add (3, "kmers_percent_weak", "%.1f", 100.0 - 100.0 * (double)nbSolids / (double)_totalKmerNb  );  }

    getInfo()->add (2, "histogram");
    getInfo()->add (3, "cutoff",            "%ld",  _histogram->get_solid_cutoff());
    getInfo()->add (3, "nb_ge_cutoff",      "%ld",  _histogram->get_nbsolids_auto());
    getInfo()->add (3, "percent_ge_cutoff", "%.1f", nbSolids > 0 ? 100.0 * (double)_histogram->get_nbsolids_auto() / (double)_bankStats.kmersNbValid : 0);
    getInfo()->add (3, "first_peak",         "%ld",  _histogram->get_first_peak());

    double N = ((double)_histogram->get_first_peak() * _bankStats.getSeqMean()) / (_bankStats.getSeqMean() - _kmerSize + 1);
    if (N > 0)  {  getInfo()->add (3, "genome_size_estimate", "%.0f",  (double)_bankStats.sequencesTotalLength / N);  }

    size_t smallestPartition = ~0;
    size_t biggestPartition  = 0;
    for (size_t i=0; i<_solidCounts->size(); i++)
    {
        size_t currentNb = (*_solidCounts)[i].getNbItems();
        smallestPartition = std::min (smallestPartition, currentNb);
        biggestPartition  = std::max (biggestPartition,  currentNb);
    }
    getInfo()->add (2, "partitions");
    getInfo()->add (3, "nb_partitions", "%ld", _solidCounts->size());
    getInfo()->add (3, "nb_items",      "%ld", _solidCounts->getNbItems());
    getInfo()->add (3, "part_biggest",  "%ld", biggestPartition);
    getInfo()->add (3, "part_smallest", "%ld", smallestPartition);
    if (_solidCounts->size())
    {
        getInfo()->add (3, "part_mean",     "%.1f", (double)nbSolids / (double)_solidCounts->size());
    }
    getInfo()->add (3, "cmd",           "hash:%d vector:%d", _partCmdTypes.first, _partCmdTypes.second);

    _fillTimeInfo /= getDispatcher()->getExecutionUnitsNumber();
    getInfo()->add (2, _fillTimeInfo.getProperties("fillsolid_time"));

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

// estimated the number of distinct kmers in a dataset
// wrapper around a Linear Counter. Adapted from Kmergenie code.
// why not a hyperloglog? it seems that the transition from the 32bit-hash initial implementation to 64bits, and supporting billions of elements, is nontrivial, so i didn't bother
// probably deserves to be in its own file
template<size_t span>
class EstimateNbDistinctKmers 
{
public:

    /** Shortcut. */
    typedef typename SortingCountAlgorithm<span>::Type  Type;
    typedef typename SortingCountAlgorithm<span>::Model Model;
    typedef typename Model::Kmer                        KmerType;

    /** */
    void estimate()
    {
        nb_distinct_kmers =(unsigned long)( (float)(linearCounter->count( )) * ((float)nbKmersTotal / (float)nbProcessedKmers)); // dubious linear extrapolation, that's all I got

        abs_error = abs((long)(nb_distinct_kmers-previous_nb_distinct_kmers));

        previous_nb_distinct_kmers = nb_distinct_kmers;
    }

    /** */
    void operator() (Sequence& sequence)
    {
        /** We build the kmers from the current sequence. */
        if (model.build (sequence.getData(), kmers) == false)  {  throw "reached EOF"; return; }

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            linearCounter->add((kmers[i].value()));
            

            // heuristics to stop early, i found that it's inaccurate with low coverage (e.g. on dsk/test/FiftyK.fastq)
            /* 
            if (nbProcessedReads % eval_every_N_reads == 0 )
            {
                
                // let's see if the estimation converges..
                // the following stopping condition will grossly over-estimate the number of distinct kmers
                // but I expect the correct result to be in the same order of magnitude
                // and better to overestimate than underestimate (for both dsk and kmergenie)
                   estimate(); 
                   bool debug = true;
                   if (debug)
                       printf("linear estimator at %ld kmers, number of distinct kmers estimated now: %ld, abs error: %ld\n",nbProcessedKmers, nb_distinct_kmers, abs_error);
                   if (abs_error < previous_nb_distinct_kmers/20) // 5% error
                   {
                       throw "LinearCounter converged"; // well, "converged" is a big word
                       return;
                   }
                   if (!linearCounter->is_accurate())
                   {
                   printf("LinearCounter is inaccurate";
                   return;
                   }

            }*/

        }
        nbProcessedKmers += kmers.size();
        nbProcessedReads++;
        //if (nbProcessedReads % 100000 == 0) printf("nb: %ld\n",nbProcessedReads);

        // disabled progress 
        //if (nbCurProgressKmers > 500000)   {  _progress.inc (nbCurProgressKmers);  nbCurProgressKmers = 0;  }
    }

    EstimateNbDistinctKmers (Model& model, u_int32_t max_memory, unsigned long nb_kmers_total, IteratorListener* progress)
        : model(model),  eval_every_N_reads(10000000),   nbKmersTotal(nb_kmers_total), 
        nbProcessedKmers(0), nbCurProgressKmers(0), previous_nb_distinct_kmers(0), nbProcessedReads(0), abs_error(0)
        //, _progress  (progress,System::thread().newSynchronizer())  
    {
        unsigned long size_linearCounter; // (in bits)
        /* let's set it to just use half of all memory available at most, ok? this isn't very robust for huge dataset, so to be tested*/
        /* if it's a tiny dataset, let's set it to total number of kmers */
        size_linearCounter = std::min( nb_kmers_total, (unsigned long) (max_memory*8*1024*1024/2) );  
        linearCounter =  new LinearCounter<span>(size_linearCounter);
    }

    unsigned long getEstimation()
    {
        estimate();
        // soo.. if it's not accurate, let's assume we have a hugeload of kmers, and let's be safe, we return the total number of kmers
        if (!linearCounter->is_accurate())
        {   
            cout << "Warning: linear counter was not accurate, returning worst-case estimation of number of distinct kmers";
            return nbKmersTotal;
        }
        return nb_distinct_kmers;
    }

private:

    /** Local resources. */
    Model&    model;
    unsigned long nbProcessedReads, nbProcessedKmers;
    unsigned long nbCurProgressKmers;
    unsigned long nbKmersTotal;
    unsigned long abs_error;
    vector<KmerType> kmers;
    LinearCounter<span> *linearCounter;
    int eval_every_N_reads;
    unsigned long previous_nb_distinct_kmers, nb_distinct_kmers;
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
void SortingCountAlgorithm<span>::configure (IBank* bank)
{
    float load_factor = 0.7;

   /** We get some information about the bank. */

    /** By default, we want to have mmers of size 8. However (for unit tests for instance),
     * we may need to have kmer sizes less than 8; in such a case, we set by convention m=k-1. */
    if (_minim_size == 0)
        _minim_size = 8;

    _minim_size = std::min ((int)_kmerSize-1, (int)_minim_size);

    // optimism == 0 mean that we guarantee worst case the memory usage,
    // any value above assumes that, on average, any distinct k-mer will be seen 'optimism+1' times
    int optimism = 0; // 0: guarantees to always work; above 0: risky
 
    /** We get some information about the bank. */
    bank->estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));
    u_int64_t bankSize = _estimateSeqTotalSize / MBYTE;

    _volume =  kmersNb * sizeof(Type) / MBYTE;  // in MBytes

    if (_volume == 0)   { _volume = 1; }    // tiny files fix

    u_int64_t volume_minim = _volume * 0.5 *1.2  ; //0.5 for using kxmers   1.2 if bad repartition of minimizers ( todo sampling to assert ram usage)

    if (volume_minim == 0)   { volume_minim = 1; }    // tiny files fix

    if (_max_disk_space == 0)  { _max_disk_space = std::min (available_space/2, 3*bankSize);  }  // used to be just bankSize until Oct 2014, changed that to 3x
    if (_max_disk_space == 0)  { _max_disk_space = 10000; }

    if (_max_memory == 0)  {  _max_memory = System::info().getMemoryProject(); }
    if (_max_memory == 0)  {  _max_memory = 1000; }
   
    /* make sure to not use more mem than system, when max_memory has default value (useful for docker images) */
    if (_max_memory == 2000)  {  
        unsigned long system_mem = System::info().getMemoryPhysicalTotal() / MBYTE;
        if (_max_memory > (system_mem * 2) / 3)
        {
            _max_memory = (system_mem * 2) / 3;
            cout << "Warning: default memory usage (2000 MB) is close or above system max, setting memory to: " << _max_memory << " MB" << endl;
        }
    }

    assert (_max_disk_space > 0);
    
    //_nb_passes = ( (_volume/3) / _max_disk_space ) + 1; //minim, approx volume /3
    _nb_passes = 1; //do not constrain nb passes on disk space anymore (anyway with minim, not very big)
    //increase it only if ram issue

    //printf("_volume  %lli volume_minim %lli _max_disk_space %lli  _nb_passes init %i  \n", _volume,volume_minim,_max_disk_space,_nb_passes);
    size_t max_open_files = System::file().getMaxFilesNumber() / 2;
    u_int64_t volume_per_pass;
    float est_volume_distinct_ratio; 

    /* disabled by default; this was an experiment */
    if (_flagEstimateNbDistinctKmers)
    {
        /* we estimate the volume of distinct kmers vs total number of kmers.
         * we store it in the variable "est_volume_distinct_ratio"
         * to compute it, we need a linear counter, let's call it now */

        TIME_INFO (getTimeInfo(), "estimate_distinct_kmers");
        Iterator<Sequence>* itSeq = _bank->iterator();
        LOCAL (itSeq);

        //_progress->setMessage (progressFormat0); // not touching progress here anymore
        Model model (_kmerSize, _minim_size);
        EstimateNbDistinctKmers<span> estimate_nb_distinct_kmers_function(model, _max_memory, kmersNb, _progress);

        /** We launch the iteration of the sequences iterator with the created functors. */
        try {
            itSeq->iterate (estimate_nb_distinct_kmers_function);
        }
        catch (const char* except)
        {

        }
        _estimatedDistinctKmerNb = estimate_nb_distinct_kmers_function.getEstimation();
        est_volume_distinct_ratio = (float) _estimatedDistinctKmerNb / (float)kmersNb;
        //est_volume_distinct_ratio = 1; // for debug
        /* est_volume_distinct_ratio == 1 mean that we guarantee worst case the memory usage,
           the value mean that, on average, a k-mer will be seen 'est_volume_distinct_ratio' times */
        // if wrongly estimated, the error 'OAHash: max rehashes..' can happen
        printf ("LinearCounter done, estimated %ld number of distinct kmers, ratio to total number of kmers: %.2f\n", (long)_estimatedDistinctKmerNb, est_volume_distinct_ratio);
    }

    do  {

        assert (_nb_passes > 0);
        volume_per_pass = volume_minim / _nb_passes;

        assert (_max_memory > 0);
        //printf("volume_per_pass %lli  _nbCores %zu _max_memory %i \n",volume_per_pass, _nbCores,_max_memory);

        if (_partitionType == 1) // adjust partition size for hash table
        {
            _nb_partitions = (u_int32_t) ceil((float) _nb_partitions / load_factor);
            _nb_partitions = ((_nb_partitions * OAHash<Type>::size_entry()) + sizeof(Type)-1) / sizeof(Type); // also adjust for hash overhead
            if (_flagEstimateNbDistinctKmers)
            {
                // use our estimation of number of distinct kmers to refine number of partitions
                // it's essentially a way to set optimism optimally
                // i'm not enabling it because computing it is slow, and reward was too small
                _nb_partitions = std::max ((u_int32_t) ceil( (float) _nb_partitions *  est_volume_distinct_ratio  * 1.3 ), (u_int32_t)1);  // 1.3 is for security
            }
            else
            {
                _nb_partitions = std::max ((_nb_partitions/(optimism+1)), (u_int32_t)1);
            }
        }
        else
        {
            // _nb_partitions  = ( (volume_per_pass*_nbCores) / _max_memory ) + 1;
            _nb_partitions  = ( ( volume_per_pass* _nb_partitions_in_parallel) / _max_memory ) + 1;

            //printf("nb passes  %i  (nb part %i / %zu)\n",_nb_passes,_nb_partitions,max_open_files);
            //_nb_partitions = max_open_files; break;
        }
        if (_nb_partitions >= max_open_files &&  _nb_partitions_in_parallel >1)   { _nb_partitions_in_parallel  = _nb_partitions_in_parallel /2;  }
        else if (_nb_partitions >= max_open_files && _nb_partitions_in_parallel == 1)   { _nb_passes++;  }
        else                                    { break;         }

        //printf("update nb passes  %i  (nb part %i / %zu)\n",_nb_passes,_nb_partitions,max_open_files);
    } while (1);

    if (_nb_partitions < 50 &&  (max_open_files - _nb_partitions  > 30) ) _nb_partitions += 30; //some more does not hurt
    
    //round nb parti to upper multiple of _nb_partitions_in_parallel if possible
    int  incpart = _nb_partitions_in_parallel - _nb_partitions % _nb_partitions_in_parallel;
    incpart = incpart % _nb_partitions_in_parallel;
    if((max_open_files - _nb_partitions  > incpart)) _nb_partitions+= incpart ;

    //_nb_partitions_in_parallel = 1 ;

    //then put _nbCores_per_partition

    _nbCores_per_partition =  _nbCores / _nb_partitions_in_parallel ; //how to best distrib available cores ?

    // with this formula we'll sometimes use less than _nbCores (maybe enforce _nb_partitions_in_parallel is power of two ?)
    DEBUG (("_nbCores %zu  _nb_partitions_in_parallel %zu   _nbCores_per_partition %zu  nb part %u nb passes %i\n",
        _nbCores,_nb_partitions_in_parallel,_nbCores_per_partition,_nb_partitions,_nb_passes
    ));

    assert(_nbCores_per_partition > 0);

    /** Now, we can define the output solid partition. */
    setSolidCounts (& (*_storage)("dsk").getPartition<Count> ("solid", _nb_partitions));

    /** We gather some statistics. */
    getInfo()->add (1, "config");
    getInfo()->add (2, "kmer_size",         "%ld", _kmerSize);
    getInfo()->add (2, "mini_size",         "%ld", _minim_size);
    getInfo()->add (2, "solidity_kind",      "%s", toString(_solidityKind).c_str());
    getInfo()->add (2, "abundance_min",     "%ld", _abundance.first);
    getInfo()->add (2, "abundance_max",     "%ld", _abundance.second);
    getInfo()->add (2, "available_space",   "%ld", available_space);
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
    getInfo()->add (2, "minimizer_type",    "%s",  (_minimizerType == 0) ? "lexicographic (kmc2 heuristic)" : "frequency");
    getInfo()->add (2, "repartition_type",    "%s",  (_repartitionType == 0) ? "unordered" : "ordered");
    if  (_flagEstimateNbDistinctKmers)
    {
        getInfo()->add (2, "estimated_nb_distinct_kmers",     "%ld", _estimatedDistinctKmerNb);
        getInfo()->add (2, "est_volume_distinct_ratio",    "%f",  est_volume_distinct_ratio);
    }
    getInfo()->add (2, "nb_cores_per_partition",     "%d",  _nbCores_per_partition);
    getInfo()->add (2, "nb_partitions_in_parallel",  "%d",  _nb_partitions_in_parallel);
}

/********************************************************************************/

/* This functor class takes a Sequence as input and split it into super kmers.
 * Each time a new superkmer is found, the 'processSuperkmer' method is called.
 *
 * NOTE : 'processSuperkmer' is pure virtual and so Sequence2SuperKmer class has to be
 * inherited for providing actual implementation of 'processSuperkmer'.
 *
 * NOTE : 'processSuperkmer' is virtual and called many times, which implies a
 * time overhead. Preliminary tests seem to show that this overhead is acceptable;
 * otherwise a template based approach could be used instead (ie Sequence2SuperKmer is
 * templated by a functor that does the job done in 'processSuperkmer').
 */
template<size_t span>
class Sequence2SuperKmer
{
public:
    /** Shortcut. */
    typedef typename SortingCountAlgorithm<span>::Type            Type;
    typedef typename SortingCountAlgorithm<span>::ModelCanonical  ModelCanonical;
    typedef typename SortingCountAlgorithm<span>::Model           Model;
    typedef typename Model::Kmer                                  KmerType;
    typedef typename Kmer<span>::SuperKmer                        SuperKmer;

    void operator() (Sequence& sequence)
    {
        /** We update statistics about the bank. */
        _bankStatsLocal.update (sequence);

        /** We first check whether we got kmers from the sequence or not. */
        if (_model.build (sequence.getData(), _kmers) == false)  { return; }

        int maxs = (Type::getSize() - 8 )/2 ;  // 8 is because  8 bit used for size of superkmers, not mini size

        /** We create a superkmer object. */
        SuperKmer superKmer (_kmersize, _miniSize, _kmers);

        /** We loop over the kmers of the sequence. */
        for (size_t i=0; i<_kmers.size(); i++)
        {
            if (_kmers[i].isValid() == false)
            {
                // on invalid kmer : output previous superk utput prev
                processSuperkmer (superKmer);

                superKmer.minimizer = DEFAULT_MINIMIZER;  //marking will have to restart 'from new'
                superKmer.range     = make_pair(i+1,i+1);

                _bankStatsLocal.kmersNbInvalid ++;

                continue;
            }

            _bankStatsLocal.kmersNbValid ++;

            /** We get the value of the current minimizer. */
            u_int64_t h = _kmers[i].minimizer().value().getVal();

            /** We have to set minimizer value if not defined. */
            if (superKmer.isValid() == false)  {  superKmer.minimizer = h;  }

            /** If the current super kmer is finished (or max size reached), we dump it. */
            if (h != superKmer.minimizer || superKmer.size() >= maxs)
            {
                processSuperkmer (superKmer);
                superKmer.range = make_pair(i,i);
            }

            /** We update the superkmer properties (minimizer value and kmers range). */
            superKmer.minimizer    = h;
            superKmer.range.second = i;
        }

        //output last superK
        processSuperkmer (superKmer);

        if (_nbWrittenKmers > 500000)   {  _progress->inc (_nbWrittenKmers);  _nbWrittenKmers = 0;  }
    }

    /** Constructor. */
    Sequence2SuperKmer (
        Model&            model,
        size_t            nbPasses,
        size_t            currentPass,
        size_t            nbPartitions,
        IteratorListener* progress,
        BankStats&        bankStats
    )
    : _model(model), _nbPass(nbPasses), _pass(currentPass), _nbPartitions(nbPartitions),
      _nbWrittenKmers(0), _progress (progress), _nbSuperKmers(0),
      _bankStatsGlobal(bankStats)
    {
        /** Shortcuts. */
        _kmersize = model.getKmerSize();
        _miniSize = model.getMmersModel().getKmerSize();
    }

    /** Destructor (virtual). */
    virtual ~Sequence2SuperKmer ()  {  _bankStatsGlobal += _bankStatsLocal;  }

protected:

    Model&           _model;
    size_t           _pass;
    size_t           _nbPass;
    size_t           _nbPartitions;
    vector<KmerType> _kmers;
    size_t           _kmersize;
    size_t           _miniSize;
    IteratorListener* _progress;
    size_t           _nbWrittenKmers;
    size_t           _nbSuperKmers;
    BankStats&       _bankStatsGlobal;
    BankStats        _bankStatsLocal;

    /** Primitive of the template method operator() */
    virtual void processSuperkmer (SuperKmer& superKmer) { _nbSuperKmers++; }
};

/********************************************************************************/
/* This functor class takes a Sequence as input, splits it into super kmers and
 * get information about the distribution of minimizers.
 */
template<size_t span>
class SampleRepart  : public Sequence2SuperKmer<span>
{
    //ie ce sera posible d avoir plus dinfo , estim ram max par exemple ?
public:

    /** Shortcut. */
    typedef typename Sequence2SuperKmer<span>::Type             Type;
    typedef typename Sequence2SuperKmer<span>::ModelCanonical   ModelCanonical;
    typedef typename Sequence2SuperKmer<span>::Model            Model;
    typedef typename Model::Kmer                           KmerType;
    typedef typename Kmer<span>::SuperKmer                 SuperKmer;

    /** */
    void processSuperkmer (SuperKmer& superKmer)
    {
        DEBUG (("SampleRepart: should count superk %i \n", superKmer.size()));

        if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid() ) //check if falls into pass
        {
            bool prev_which = superKmer[0].which();
            size_t kx_size = 0;
                    
            /** Shortcut. */
            size_t superKmerLen = superKmer.size();
            
            /** We increase superkmer counter the current minimizer. */
            _local_pInfo.incSuperKmer_per_minimBin (superKmer.minimizer, superKmerLen);

            /** We loop over the kmer of the superkmer (except the first one).
             *  We update the pInfo each time we find a kxmer in the superkmer. */
            for (size_t ii=1 ; ii < superKmerLen; ii++)
            {
                /** A kxmer is defined by having successive canonical kmers. Here, we just care that
                 * successive kmer values are on the same strand. */
                if (superKmer[ii].which() != prev_which || kx_size >= _kx) // kxmer_size = 1 //cost should diminish with larger kxmer
                {
                    /** We increase the number of kxmer found for the current minimizer. */
                    _local_pInfo.incKxmer_per_minimBin (superKmer.minimizer);
                    kx_size = 0;
                }
                else
                {
                    kx_size++;
                }

                prev_which = superKmer[ii].which() ;
            }

            /** We add the pending kxmer to the bin. */
            _local_pInfo.incKxmer_per_minimBin (superKmer.minimizer);
        }
    }

    /** Constructor. */
    SampleRepart (
        Model&            model,
        size_t            nbPasses,
        size_t            currentPass,
        size_t            nbPartitions,
        IteratorListener* progress,
        BankStats&        bankStats,
        PartiInfo<5>&     pInfo
    )
    :   Sequence2SuperKmer<span> (model, nbPasses, currentPass, nbPartitions, progress, bankStats),
        _kx(4), _extern_pInfo(pInfo), _local_pInfo(nbPartitions,model.getMmersModel().getKmerSize())
    {
    }

    /** Destructor. */
    ~SampleRepart ()
    {
        //add to global parti_info
        _extern_pInfo += _local_pInfo;
    }

private:
    size_t        _kx;
    PartiInfo<5>& _extern_pInfo;
    PartiInfo<5>  _local_pInfo;
};


template<size_t span>
class MmersFrequency
{
public:
    /** Shortcut. */
    typedef typename SortingCountAlgorithm<span>::Type            Type;
    typedef typename SortingCountAlgorithm<span>::ModelCanonical  ModelCanonical;
    typedef typename SortingCountAlgorithm<span>::Model           Model;
    typedef typename Model::Kmer                                  KmerType;

	typedef typename Kmer<span>::ModelDirect     ModelDirect;
	typedef typename ModelDirect::Kmer     KmerTypeDirect;

    void operator() (Sequence& sequence)
    {
        /** We first check whether we got mmers from the sequence or not. */
        if (_minimodel->build (sequence.getData(), _mmers) == false)  { return; }

        /** We loop over the mmers of the sequence. */
        for (size_t i=0; i<_mmers.size(); i++)
        {
            if (_mmers[i].isValid() == false)
                continue;

            /** increment m-mer count */
            _m_mer_counts[_mmers[i].value().getVal()] ++;
        }

        if (_nbProcessedMmers > 500000)   {  _progress.inc (_nbProcessedMmers);  _nbProcessedMmers = 0;  }
    }

    /** Constructor. */
    MmersFrequency (
        int mmerSize,
        IteratorListener* progress,
        uint32_t*         m_mer_counts
    )
    : 
      _nbProcessedMmers(0), _progress (progress,System::thread().newSynchronizer()),
      _m_mer_counts(m_mer_counts)
    {
        _minimodel = new ModelDirect(mmerSize); // FIXME: should it be ModelCanonical??
        u_int64_t nbminim = (uint64_t)pow(4.0,mmerSize);

        for (int i = 0; i < nbminim; i++)
            _m_mer_counts[i] = 0;
    }

protected:

    ModelDirect*           _minimodel;
    vector<KmerTypeDirect> _mmers;
    size_t           _mmersize;
    ProgressSynchro  _progress;
    uint32_t*        _m_mer_counts;
    size_t           _nbProcessedMmers;
};


/********************************************************************************/
/* This functor class takes a Sequence as input, splits it into super kmers and
 * serialize them into partitions.
 *
 * A superkmer is dumped into a partition 'p' chosen by some hash code of the minimizer
 * of the superkmer. Such a hash code can be computed in several way; actually, we use
 * a lookup table that has computed the minimizers distribution on a subset of the
 * processed bank.
 */
template<size_t span>
class FillPartitions : public Sequence2SuperKmer<span>
{
public:
    /** Shortcut. */
    typedef typename Sequence2SuperKmer<span>::Type            Type;
    typedef typename Sequence2SuperKmer<span>::ModelCanonical  ModelCanonical;
    typedef typename Sequence2SuperKmer<span>::Model           Model;
    typedef typename Model::Kmer                          KmerType;
    typedef typename Kmer<span>::SuperKmer                SuperKmer;

    /** */
    void processSuperkmer (SuperKmer& superKmer)
    {
        if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid()) //check if falls into pass
        {
            /** We get the hash code for the current miminizer.
             * => this will give us the partition where to dump the superkmer. */
            size_t p = this->_repartition (superKmer.minimizer);

            /** We save the superkmer into the right partition. */
            superKmer.save (this->_partition[p]);

            /*********************************************/
            /** Now, we compute statistics about kxmers. */
            /*********************************************/

            Type radix, radix_kxmer_forward ,radix_kxmer ;
            bool prev_which = superKmer[0].which();
            size_t kx_size =0;

            radix_kxmer_forward = getHeavyWeight (superKmer[0].value());

            for (size_t ii=1 ; ii < superKmer.size(); ii++)
            {
                //compute here stats on  kx mer
                //tant que tai <= xmer et which kmer[ii] == which kmer [ii-1] --> cest un kxmer
                //do the same in sampling : gives ram estimation
                if (superKmer[ii].which() != prev_which || kx_size >= this->_kx) // kxmer_size = 1 //cost should diminish with larger kxmer
                {
                    //output kxmer size kx_size,radix_kxmer
                    //kx mer is composed of _superKp[ii-1] _superKp[ii-2] .. _superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
                    if(prev_which)
                    {
                        radix_kxmer = radix_kxmer_forward;
                    }
                    else // si revcomp, le radix du kxmer est le debut du dernier kmer
                    {
                        radix_kxmer = getHeavyWeight (superKmer[ii-1].value());
                    }

                    this->_local_pInfo.incKmer_and_rad (p, radix_kxmer.getVal(), kx_size); //nb of superkmer per x per parti per radix

                    radix_kxmer_forward =  getHeavyWeight (superKmer[ii].value());
                    kx_size =0;
                }
                else
                {
                    kx_size++;
                }

                prev_which = superKmer[ii].which() ;
            }

            //record last kx mer
            if(prev_which)
            {
                radix_kxmer = radix_kxmer_forward;
            }
            else // si revcomp, le radix du kxmer est le debut du dernier kmer
            {
                radix_kxmer =  getHeavyWeight (superKmer[superKmer.size()-1].value());
            }

            this->_local_pInfo.incKmer_and_rad(p, radix_kxmer.getVal(),kx_size );

            /** We update progression information. */
            this->_nbWrittenKmers += superKmer.size();
        }
    }

    /** Constructor. */
    FillPartitions (
        Model&             model,
        size_t             nbPasses,
        size_t             currentPass,
        size_t             nbPartitions,
        IteratorListener*  progress,
        BankStats&         bankStats,
        Partition<Type>*   partition,
        Repartitor&        repartition,
        PartiInfo<5>&      pInfo
    )
    :   Sequence2SuperKmer<span> (model, nbPasses, currentPass, nbPartitions, progress, bankStats),
        _kx(4),
        _extern_pInfo(pInfo) , _local_pInfo(nbPartitions,model.getMmersModel().getKmerSize()),
        _repartition (repartition), _partition (*partition,1<<12,0)
    {
        _mask_radix = (int64_t) 255 ;
        _mask_radix = _mask_radix << ((this->_kmersize - 4)*2); //get first 4 nt  of the kmers (heavy weight)
    }

    /** Destructor. */
    ~FillPartitions ()
    {
        //add to global parti_info
        _extern_pInfo += _local_pInfo;
    }

private:

    size_t        _kx;
    PartiInfo<5>& _extern_pInfo;
    PartiInfo<5>  _local_pInfo;
    Type          _mask_radix;
    Repartitor&   _repartition;

    /** Shared resources (must support concurrent accesses). */
#ifdef PROTO_COMP
    PartitionCacheSorted<Type> _partition;
#else
    PartitionCache<Type> _partition;
#endif

    Type getHeavyWeight (const Type& kmer) const  {  return (kmer & this->_mask_radix) >> ((this->_kmersize - 4)*2);  }
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
void SortingCountAlgorithm<span>::fillPartitions (size_t pass, Iterator<Sequence>* itSeq, PartiInfo<5>& pInfo)
{
    TIME_INFO (getTimeInfo(), "fill_partitions");

    DEBUG (("SortingCountAlgorithm<span>::fillPartitions  _kmerSize=%d _minim_size=%d \n", _kmerSize, _minim_size));

    /** We delete the previous partitions storage. */
    if (_partitionsStorage)  { _partitionsStorage->remove (); }

    /** We build the temporary storage name from the output storage name. */
    string tmpStorageName = _storage->getName() + string("_partitions");

    /** We create the partition files for the current pass. */
#ifdef PROTO_COMP
    setPartitionsStorage (StorageFactory(STORAGE_COMPRESSED_FILE).create ("partitions", true, false));
    //setPartitionsStorage (StorageFactory(STORAGE_GZFILE).create ("partitions", true, false));
#else
    setPartitionsStorage (StorageFactory(STORAGE_FILE).create (tmpStorageName, true, false));
#endif
    
    setPartitions        (0); // close the partitions first, otherwise new files are opened before  closing parti from previous pass

    setPartitions        ( & (*_partitionsStorage)().getPartition<Type> ("parts", _nb_partitions));

    /** We update the message of the progress bar. */
    _progress->setMessage (Stringify::format(progressFormat1, _current_pass+1, _nb_passes));

    u_int64_t nbseq_sample = std::max ( u_int64_t (_estimateSeqNb * 0.05) ,u_int64_t( 1000000ULL) ) ;

    DEBUG (("SortingCountAlgorithm<span>::fillPartitions : nb seq for sample :  %llu \n ",nbseq_sample));

    /* now is a good time to switch to frequency-based minimizers if required:
      because right after we'll start using minimizers to compute the distribution 
      of superkmers in bins */
    uint32_t *freq_order = NULL;
    std::vector<std::pair<int, int> > counts;
    if (_minimizerType == 1)
    {
        u_int64_t rg = ((u_int64_t)1 << (2*_minim_size));
        //cout << "\nAllocating " << ((rg*sizeof(uint32_t))/1024) << " KB for " << _minim_size <<"-mers frequency counting (" << rg << " elements total)" << endl;
        uint32_t *m_mer_counts = new uint32_t[rg];
        Model model( _kmerSize,_minim_size);

        // can we reuse the it_sample variable above?
        Iterator<Sequence>* it_sample = createIterator (
                new TruncateIterator<Sequence> (*itSeq, nbseq_sample),
                nbseq_sample,
                "Approximating frequencies of minimizers" 
                );
        LOCAL (it_sample);

        /** We compute an estimation of minimizers frequencies from a part of the bank. */
        // actually.. let's try with the whole thing (itSeq instead of it_sample)
        getDispatcher()->iterate (it_sample,  MmersFrequency<span> (
            _minim_size, _progress, m_mer_counts)
        );
       
        // single threaded, for debugging
        /*MmersFrequency<span> mmersfrequency(model, _progress, bstatsDummy, m_mer_counts);
        it_sample->iterate(mmersfrequency);*/

        /* sort frequencies */
        for (int i(0); i < rg; i++)
        {
            if (m_mer_counts[i] > 0)
                counts.push_back(make_pair(m_mer_counts[i],i));
        }
        delete[] m_mer_counts;

        sort(counts.begin(),counts.end());

        /* assign frequency to minimizers */
        freq_order = new uint32_t[rg];

        for (int i = 0; i < rg ; i++)
            freq_order[i] = rg; // set everything not seen to highest value (not a minimizer)

        for (unsigned int i = 0; i < counts.size(); i++)
        {
            freq_order[counts[i].second] = i;
        }

        // small but necessary trick: the largest minimizer has to have largest rank, as it's used as the default "largest" value 
        freq_order[rg-1] = rg-1;

        model.setMinimizersFrequency(freq_order);
   
        // save this function 
        tools::storage::impl::Storage::ostream os (getStorageGroup(), "minimFrequency");
        os.write ((const char*)freq_order,    sizeof(uint32_t) * rg);
        os.flush();
    }

    /** We create a kmer model; using the frequency order if we're in that mode */
    Model model( _kmerSize,_minim_size, typename kmer::impl::Kmer<span>::ComparatorMinimizerFrequency(), freq_order);

    int mmsize = model.getMmersModel().getKmerSize();

    PartiInfo<5> sample_info (_nb_partitions,mmsize);

    /** We create an iterator over a truncated part of the input bank. */
    Iterator<Sequence>* it_sample = createIterator (
        new TruncateIterator<Sequence> (*itSeq, nbseq_sample),
        nbseq_sample,
        progressFormat3
    );
    LOCAL (it_sample);

    BankStats bstatsDummy;

    /** We compute a distribution of Superkmers from a part of the bank. */
    getDispatcher()->iterate (it_sample,  SampleRepart<span> (
        model, _nb_passes, pass, _nb_partitions, _progress, bstatsDummy, sample_info)
    );

    /** We compute the distribution of the minimizers. As a result, we will have a hash function
     * that gives a hash code for a minimizer value. */
    Repartitor repartitor (_partitions->size(), mmsize);
    if (_minimizerType == 1)
        repartitor.justGroup (sample_info, counts);
    else
    {
        repartitor.computeDistrib (sample_info);
        if (_repartitionType == 1)
        {
            repartitor.justGroupLexi (sample_info); // For bcalm, i need the minimizers to remain in order. so using this suboptimal but okay repartition
        }
        
    }

    /** We save the distribution (may be useful for debloom for instance). */
    repartitor.save (getStorageGroup());

	/** We have to reinit the progress instance since it may have been used by SampleRepart before. */
    _progress->init();

    /** We may have several input banks instead of a single one. */
    std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();

    /** We first reset the vector holding the kmers number for each partition and for each bank.
     * It can be seen as the following matrix:
     *
     *           part0  part1  part2 ... partJ
     *   bank0    xxx    xxx    xxx       xxx
     *   bank1    xxx    xxx    xxx       xxx
     *    ...
     *   bankI    xxx    xxx    xxx       xxx
     *
     *   Here xxx is the number of items found for the bank I in the partition J
     */
    _nbKmersPerPartitionPerBank.clear();

    /** We launch the iteration of the sequences iterator with the created functors. */
    for (size_t i=0; i<itBanks.size(); i++)
    {
        /** We fill the partitions. */
        getDispatcher()->iterate (itBanks[i], FillPartitions<span> (
            model, _nb_passes, pass, _nb_partitions, _progress, _bankStats, _partitions, repartitor, pInfo
        ));

        /** We flush the partitions in order to be sure to have the exact number of items per partition. */
        _partitions->flush();

        /** We get a snapshot of items number in each partition. */
        vector<size_t> nbItems;
        for (size_t p=0; p<_nb_partitions; p++)
        {
            nbItems.push_back ((*_partitions)[p].getNbItems());
        }

        /** We add the current number of kmers in each partition for the reached ith bank. */
        _nbKmersPerPartitionPerBank.push_back (nbItems);
    }
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
std::vector<size_t> SortingCountAlgorithm<span>::getNbCoresList (PartiInfo<5>& pInfo)
{
    std::vector<size_t> result;

    for (size_t p=0; p<_nb_partitions; )
    {
        u_int64_t ram_total = 0;
        size_t i=0;
        for (i=0; i< _nb_partitions_in_parallel && p<_nb_partitions
            && (ram_total ==0  || ((ram_total+(pInfo.getNbSuperKmer(p)*getSizeofPerItem()))  <= _max_memory*MBYTE)) ; i++, p++)
        {
            ram_total += pInfo.getNbSuperKmer(p)*getSizeofPerItem();
        }

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
void SortingCountAlgorithm<span>::fillSolidKmers (PartiInfo<5>& pInfo)
{
    TIME_INFO (getTimeInfo(), "fill_solid_kmers");

    DEBUG (("SortingCountAlgorithm<span>::fillSolidKmers\n"));

    /** We update the message of the progress bar. */
    _progress->setMessage (Stringify::format (progressFormat2, _current_pass+1, _nb_passes));

    /** We retrieve the list of cores number for dispatching N partitions in N threads.
     *  We need to know these numbers for allocating the N maps according to the maximum allowed memory.
     */
    vector<size_t> coreList = getNbCoresList(pInfo); //uses _nb_partitions_in_parallel

    /** We need a memory allocator. */
    MemAllocator pool;

    size_t p = 0;
    for (size_t i=0; i<coreList.size(); i++)
    {
        vector<ICommand*> cmds;

        size_t currentNbCores = coreList[i];
        assert (currentNbCores > 0);

        /** We correct the number of memory per map according to the max allowed memory.
         * Note that _max_memory has initially been divided by the user provided cores number. */
        u_int64_t mem = (_max_memory*MBYTE)/currentNbCores;

        /** We need to cache the solid kmers partitions.
         *  NOTE : it is important to save solid kmers by big chunks (ie cache size) in each partition.
         *  Indeed, if we directly iterate the solid kmers through a Partition::iterator() object,
         *  one partition is iterated after another one, which doesn't reflect the way they are in filesystem,
         *  (ie by chunks of solid kmers) which may lead to many moves into the global HDF5 file.
         *  One solution is to make sure that the written chunks of solid kmers are big enough: here
         *  we accept to provide at most 2% of the max memory, or chunks of 200.000 items.
         */
        size_t cacheSize = std::min ((u_int64_t)(200*1000), mem/(50*sizeof(Count)));

        DEBUG (("SortingCountAlgorithm::fillSolidKmers:  computing %zu partitions simultaneously , parti : ",currentNbCores));

        for (size_t j=0; j<currentNbCores; j++, p++)
        {
            ISynchronizer* synchro = System::thread().newSynchronizer();
            LOCAL (synchro);

            /** We get the pth collection for storing solid kmers for the current partition. */
            Bag<Count>* solidKmers = & (*_solidCounts)[p];

            DEBUG ((" %zu ", p));

            /* Get the memory taken by this partition if loaded for sorting */
            uint64_t memoryPartition = (pInfo.getNbSuperKmer(p)*getSizeofPerItem()); //in bytes
            DEBUG (("  (%llu  MB) ",memoryPartition/MBYTE));

            bool forceHashing = (_partitionType == 1);

            /** If we have several input banks, we may have to compute kmer solidity for each bank, which
             * can be currently done only with sorted vector. */
            bool forceVector  = _nbKmersPerPartitionPerBank.size() > 1;

            ICommand* cmd = 0;

            //still use hash if by vector would be too large even with single part at a time
            if ( ((memoryPartition > mem && currentNbCores==1) || forceHashing) && !forceVector)
            {
                if (pool.getCapacity() != 0)  {  pool.reserve(0);  }

                // also allow to use mem pool for oahash ? ou pas la peine
                cmd = new PartitionsByHashCommand<span>   (
                    solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb, _abundance, _progress, _fillTimeInfo,
                    pInfo,p,_nbCores_per_partition, _kmerSize,pool, cacheSize, mem
                );

                _partCmdTypes.first ++;
            }
            else
            {
                u_int64_t memoryPoolSize = _max_memory*MBYTE;

                /** In case of forcing sorted vector (multiple banks counting for instance), we may have a
                 * partition bigger than the max memory. */
                if (forceVector  &&  pInfo.getNbSuperKmer(p)*getSizeofPerItem() >= memoryPoolSize)
                {
                    static const int EXCEED_FACTOR = 2;

                    if (pInfo.getNbSuperKmer(p)*getSizeofPerItem()  < EXCEED_FACTOR*memoryPoolSize)
                    {
                        /** We accept in this case to exceed the allowed memory. */
                        memoryPoolSize = pInfo.getNbSuperKmer(p)*getSizeofPerItem();
                    }
                    else
                    {
                        /** We launch an exception. */
                        throw Exception ("memory issue: %lld required and %lld available",
                            pInfo.getNbSuperKmer(p)*getSizeofPerItem(), memoryPoolSize
                        );
                    }
                }

               //if capa pool ==0, reserve max memo , pass pool to partibyvec, will be used  for vec kmers
                if (pool.getCapacity() == 0)  {  pool.reserve (memoryPoolSize); }

                /** Recall that we got the following matrix in _nbKmersPerPartitionPerBank
                 *
                 *           part0  part1  part2 ... partJ
                 *   bank0    xxx    xxx    xxx       xxx
                 *   bank1    xxx    xxx    xxx       xxx
                 *    ...
                 *   bankI    xxx    xxx    xxx       xxx
                 *
                 *   Now, for the current partition p, we want the number of items found for each bank.
                 *
                 *              bank0   bank1   ...   bankI
                 *   offsets :   xxx     xxx           xxx
                 */
                vector<size_t> nbItemsPerBankPerPart;
                for (size_t i=0; i<_nbKmersPerPartitionPerBank.size(); i++)
                {
                    nbItemsPerBankPerPart.push_back (_nbKmersPerPartitionPerBank[i][p] - (i==0 ? 0 : _nbKmersPerPartitionPerBank[i-1][p]) );
                }

                cmd = new PartitionsByVectorCommand<span> (
                    solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb, _abundance, _progress, _fillTimeInfo,
                    pInfo,p,_nbCores_per_partition, _kmerSize, pool, cacheSize, _solidityKind, nbItemsPerBankPerPart
                );

                _partCmdTypes.second ++;
            }

            cmds.push_back (cmd);

        } /* end of for (size_t j=0; j<currentNbCores... */

        DEBUG (("\n"));

        /** We launch the commands through a dispatcher. */
        getDispatcher()->dispatchCommands (cmds, 0);

        // free internal memory of pool here
        pool.free_all();
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
