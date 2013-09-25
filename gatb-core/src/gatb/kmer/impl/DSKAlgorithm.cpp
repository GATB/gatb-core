/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <gatb/kmer/impl/DSKAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/math/NativeInt64.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

#include <math.h>
#include <algorithm>

#ifdef OMP
#include <omptl/omptl_algorithm>
#endif

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
static const char* progressFormat1 = "DSK: Pass %d/%d, Step 1: partitioning  ";
static const char* progressFormat2 = "DSK: Pass %d/%d, Step 2: counting kmers";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
DSKAlgorithm<ProductFactory,T>::DSKAlgorithm (
    Product<ProductFactory>& product,
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
    _product(product),
    _bank(0),
    _kmerSize(kmerSize), _nks(nks),
    _partitionType(partitionType), _nbCores(nbCores), _prefix(prefix),
    _progress (0),
    _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
    _max_disk_space(max_disk_space), _max_memory(max_memory), _volume(0), _nb_passes(0), _nb_partitions(0), _current_pass(0),
    _histogram (0), _histogramUri(histogramUri),
    _partitionsProduct(0), _partitions(0), _totalKmerNb(0)
{
    setBank (bank);

    /** We create the collection corresponding to the solid kmers output. */
    setSolidKmers (& product("dsk").template getCollection<Kmer<T> > ("solid"));

    /** We set the histogram instance. */
    setHistogram (new Histogram  (10000, & product("dsk").template getCollection<Histogram::Entry>("histogram") ));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
DSKAlgorithm<ProductFactory,T>::~DSKAlgorithm ()
{
    setProgress          (0);
    setBank              (0);
    setPartitionsProduct (0);
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
template<typename ProductFactory, typename T>
void DSKAlgorithm<ProductFactory,T>::execute ()
{
    /** We retrieve the actual number of cores. */
    _nbCores = getDispatcher()->getExecutionUnitsNumber();
    assert (_nbCores > 0);

    /** We set the max memory according to the number of used cores. */
    _max_memory /= _nbCores;

    /** We configure dsk by computing the number of passes and partitions we will have
     * according to the allowed disk and memory space. */
    configure (_bank);

    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = _bank->iterator();
    LOCAL (itSeq);

    /** We configure the progress bar. */
    setProgress ( createIteratorListener (2 * _volume * MBYTE / sizeof(T), "counting kmers"));
    _progress->init ();

    /*************************************************************/
    /*                       DSK MAIN LOOP                       */
    /*************************************************************/
    /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
    for (_current_pass=0; _current_pass<_nb_passes; _current_pass++)
    {
        DEBUG (("DSKAlgorithm<ProductFactory,T>::execute  pass [%ld,%d] \n", _current_pass+1, _nb_passes));

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

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "kmers_nb_valid",     "%ld", _totalKmerNb);
    getInfo()->add (2, "kmers_nb_weak",      "%ld", _totalKmerNb - _solidKmers->iterable()->getNbItems() );
    getInfo()->add (2, "kmers_nb_solid",     "%ld", _solidKmers->iterable()->getNbItems() );
    getInfo()->add (2, "kmers_percent_weak", "%.1f", 100.0 - 100.0 * (double)_solidKmers->iterable()->getNbItems() / (double)_totalKmerNb  );

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
template<typename ProductFactory, typename T>
void DSKAlgorithm<ProductFactory,T>::configure (IBank* bank)
{
    float load_factor = 0.7;

    /** We get some information about the bank. */
    bank->estimate (_estimateSeqNb, _estimateSeqTotalSize, _estimateSeqMaxSize);

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    u_int64_t kmersNb  = (_estimateSeqTotalSize - _estimateSeqNb * (_kmerSize-1));
    u_int64_t bankSize = _estimateSeqTotalSize / MBYTE;

    _volume = kmersNb * sizeof(T) / MBYTE;  // in MBytes
    if (_volume == 0)   { _volume = 1; }    // tiny files fix

    if (_max_disk_space == 0)  { _max_disk_space = std::min (available_space/2, bankSize);  }
    if (_max_disk_space == 0)  { _max_disk_space = 10000; }

    _nb_passes = ( _volume / _max_disk_space ) + 1;

    size_t max_open_files = System::file().getMaxFilesNumber() / 2;
    u_int64_t volume_per_pass;

    do  {
        volume_per_pass = _volume / _nb_passes;
        _nb_partitions  = ( volume_per_pass / _max_memory ) + 1;

        if (_partitionType == 0)
        {
            _nb_partitions = (u_int32_t) ceil((float) _nb_partitions / load_factor);
            _nb_partitions = ((_nb_partitions * OAHash<T>::size_entry()) + sizeof(T)-1) / sizeof(T); // also adjust for hash overhead
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
    getInfo()->add (2, "nb_bits_per_kmer",  "%d",  T::getSize());
    getInfo()->add (2, "nb_cores",          "%d",  getDispatcher()->getExecutionUnitsNumber());
    getInfo()->add (2, "partition_type",    "%d",  _partitionType);
}

/********************************************************************************/

template<typename ProductFactory, typename T>
class FillPartitions : public IteratorFunctor
{
public:

    void operator() (Sequence& sequence)
    {
        /** By default, we will use the provided data with a ASCII encoding. */
        Data* data = & sequence.getData();

        /** We may have to expand the binary data to integer format. */
        if (sequence.getData().getEncoding() == Data::BINARY)
        {
            size_t expandedLen = sequence.getData().size() ;

            if (expandedLen > binaryData.size())  {  binaryData.resize (expandedLen + 4);  }

            /** We convert the provided binary data into integer encoding. */
            Data::convert (sequence.getData(), binaryData);

            /** The actual data will be the binary data. */
            data = &binaryData;
        }

        /** We build the kmers from the current sequence. */
        if (model.build (*data, kmers) == false)  { return; }

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            /** We hash the current kmer. */
            T h = oahash (kmers[i]);

            /** We check whether this kmer has to be processed during the current pass. */
            if ((h % nbPass) != pass)  { continue; }

            T reduced_kmer = h / nbPass;

            /** We compute in which partition this kmer falls into. */
            size_t p = reduced_kmer % nbPartitions;

            /** We write the kmer into the bag. */
            _partition[p].insert (kmers[i]);

            nbWrittenKmers++;
        }

        if (nbWrittenKmers > 500000)   {  _progress.inc (nbWrittenKmers);  nbWrittenKmers = 0;  }
    }

    FillPartitions (Model<T>& model, size_t nbPasses, size_t currentPass, Partition<ProductFactory, T>* partition, IteratorListener* progress)
        : model(model), pass(currentPass), nbPass(nbPasses), nbPartitions(partition->size()), nbWrittenKmers(0),
          _partition(*partition,1<<12,this->newSynchro()), _progress (progress,this->newSynchro())  {}

private:

    /** Local resources. */
    Model<T>& model;
    size_t    pass;
    size_t    nbPass;
    size_t    nbPartitions;
    size_t    nbWrittenKmers;
    Data      binaryData;
    vector<T> kmers;

    /** Shared resources (must support concurrent accesses). */
    PartitionCache<ProductFactory,T> _partition;
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
template<typename ProductFactory, typename T>
void DSKAlgorithm<ProductFactory,T>::fillPartitions (size_t pass, Iterator<Sequence>* itSeq)
{
    TIME_INFO (getTimeInfo(), "fill_partitions");

    DEBUG (("DSKAlgorithm<ProductFactory,T>::fillPartitions  pass \n", pass));

    /** We create a kmer model. */
    Model<T> model (_kmerSize);

    /** We delete the previous partitions product. */
    if (_partitionsProduct)  { _partitionsProduct->remove (); }

    /** We create the partition files for the current pass. */
    setPartitionsProduct (PartitionFactory::createProduct ("partitions", false));
    setPartitions        ( & (*_partitionsProduct)().template getPartition<T> ("parts", _nb_partitions));

    /** We update the message of the progress bar. */
    _progress->setMessage (progressFormat1, _current_pass+1, _nb_passes);

    /** We launch the iteration of the sequences iterator with the created functors. */
    getDispatcher()->iterate (itSeq, FillPartitions<PartitionFactory,T> (model, _nb_passes, pass, _partitions, _progress), 15*1000);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
class PartitionsCommand : public ICommand
{
public:
    PartitionsCommand (
        DSKAlgorithm<ProductFactory,T>& algo,
        Bag<Kmer<T> >* solidKmers,
        Iterable<T>&   partition,
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
    BagCache<Kmer<T> >  _solidKmers;
    Iterable<T>&        _partition;
    HistogramCache      _histogram;
    ProgressSynchro     _progress;
    u_int64_t           _totalKmerNb;
    u_int64_t&          _totalKmerNbRef;

    void insert (const Kmer<T>& kmer)
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
template<typename ProductFactory, typename T>
class PartitionsByHashCommand : public PartitionsCommand<ProductFactory, T>
{
public:

    PartitionsByHashCommand (
        DSKAlgorithm<ProductFactory,T>& algo,
        Bag<Kmer<T> >*  solidKmers,
        Iterable<T>&    partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        u_int64_t       hashMemory
    )
        : PartitionsCommand<ProductFactory, T> (algo, solidKmers, partition, histogram, synchro, totalKmerNbRef), _hashMemory(hashMemory)  {}

    void execute ()
    {
        size_t count=0;

        /** We need a map for storing part of solid kmers. */
        OAHash<T> hash (_hashMemory);

        /** We directly fill the vector from the current partition file. */
        Iterator<T>* it = this->_partition.iterator();  LOCAL(it);

        for (it->first(); !it->isDone(); it->next())
        {
            hash.increment (it->item());

            /** Some display. */
            if (++count == 100000)  {  this->_progress.inc (count);  count=0; }
        }

        /** We loop over the solid kmers map. */
        Iterator < Abundance<T> >* itKmerAbundance = hash.iterator();
        LOCAL (itKmerAbundance);

        for (itKmerAbundance->first(); !itKmerAbundance->isDone(); itKmerAbundance->next())
        {
            /** We may add this kmer to the solid kmers bag. */
           this->insert ((Kmer<T>&) itKmerAbundance->item());
        }
    }

private:
    u_int64_t _hashMemory;
};

/********************************************************************************/
/** */
template<typename ProductFactory, typename T>
class PartitionsByVectorCommand : public PartitionsCommand<ProductFactory, T>
{
    vector<T> kmers;

public:

    PartitionsByVectorCommand (
        DSKAlgorithm<ProductFactory,T>& algo,
        Bag<Kmer<T> >*  solidKmers,
        Iterable<T>&    partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef
    )
        : PartitionsCommand<ProductFactory, T> (algo, solidKmers, partition, histogram, synchro, totalKmerNbRef)
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
        Iterator<T>* it = this->_partition.iterator();  LOCAL (it);
        size_t idx = 0;
        for (it->first(); !it->isDone(); it->next(), idx++) { kmers[idx] = it->item(); }

        /** We set the extra item to a max value, so we are sure it will sorted at the last location.
         * This trick allows to avoid extra treatment after the loop that computes the kmers abundance. */
        kmers[partitionLen] = ~0;

        /** We sort the vector. */
#ifdef OMP
        omptl::sort (kmers.begin (), kmers.end ());
#else
        std::sort (kmers.begin (), kmers.end ());
#endif

        u_int32_t abundance = 0;
        T previous_kmer = kmers.front();

        /** We loop over the sorted solid kmers. */
        for (typename vector<T>::iterator itKmers = kmers.begin(); itKmers != kmers.end(); ++itKmers)
        {
            if (*itKmers == previous_kmer)  {   abundance++;  }
            else
            {
                this->insert (Kmer<T> (previous_kmer, abundance) );

                abundance     = 1;
                previous_kmer = *itKmers;
            }
        }

        /** We update the progress bar. */
        this->_progress.inc (kmers.size());
    }
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, typename T>
std::vector<size_t> DSKAlgorithm<ProductFactory,T>::getNbCoresList ()
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
template<typename ProductFactory, typename T>
void DSKAlgorithm<ProductFactory,T>::fillSolidKmers (Bag<Kmer<T> >*  solidKmers)
{
    TIME_INFO (getTimeInfo(), "fill_solid_kmers");

    DEBUG (("DSKAlgorithm<ProductFactory,T>::fillSolidKmers\n"));

    /** We update the message of the progress bar. */
    _progress->setMessage (progressFormat2, _current_pass+1, _nb_passes);

    ISynchronizer* synchro = System::thread().newSynchronizer();

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
        u_int64_t mem = (_max_memory*MBYTE*_nbCores)/currentNbCores;

        for (size_t j=0; j<currentNbCores; j++, p++)
        {
            ICommand* cmd = 0;

            if (_partitionType == 0)
            {
                cmd = new PartitionsByHashCommand<ProductFactory, T>   (*this, solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb, mem);
            }
            else
            {
                cmd = new PartitionsByVectorCommand<ProductFactory, T> (*this, solidKmers, (*_partitions)[p], _histogram, synchro, _totalKmerNb);
            }

            cmds.push_back (cmd);
        }

        getDispatcher()->dispatchCommands (cmds, 0);
    }

    /** Some cleanup. */
    delete synchro;
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DSKAlgorithm <ProductFileFactory, NativeInt64>;
#ifdef INT128_FOUND
template class DSKAlgorithm <ProductFileFactory, NativeInt128>;
#else
template class DSKAlgorithm <ProductFileFactory, LargeInt<2> >;
#endif

template class DSKAlgorithm <ProductFileFactory, LargeInt<3> >;
template class DSKAlgorithm <ProductFileFactory, LargeInt<4> >;

/********************************************************************************/

template class DSKAlgorithm <ProductHDF5Factory, NativeInt64>;
#ifdef INT128_FOUND
template class DSKAlgorithm <ProductHDF5Factory, NativeInt128>;
#else
template class DSKAlgorithm <ProductHDF5Factory, LargeInt<2> >;
#endif

template class DSKAlgorithm <ProductHDF5Factory, LargeInt<3> >;
template class DSKAlgorithm <ProductHDF5Factory, LargeInt<4> >;


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
