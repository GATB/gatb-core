#include "bcalm_algo.hpp"

#include <libgen.h> // for basename()
#include "logging.hpp"
#include "ograph.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <tuple>

#include <gatb/tools/designpattern/impl/Command.hpp>

#include <atomic>
#include <thread>

#include "ThreadPool.h"

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/storage/impl/StorageTools.hpp>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/kmer/impl/PartiInfo.hpp>   // for repartitor 
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#define get_wtime() chrono::system_clock::now()
#ifndef diff_wtime
#define diff_wtime(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count()
#endif

//#define BINSEQ // "graph4 is not ready" according to antoine. also, initBinSeq provokes segfault at end of bcalm

#ifdef BINSEQ
#include "binSeq.h"
#define BUCKET_STR_TYPE binSeq
#define TO_BUCKET_STR(x) binSeq(x)
#define FROM_BUCKET_STR(x) (x.str())
#else
#define BUCKET_STR_TYPE string
#define TO_BUCKET_STR(x) x
#define FROM_BUCKET_STR(x) x
#endif


// timing-related variables

#define THREAD_SAFE_TIMING
#ifdef THREAD_SAFE_TIMING
typedef std::atomic<double> atomic_double;
#else
#define atomic_double_add(d1,d2) d1 += d2;
typedef double atomic_double;
#endif


using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;






/*
 * some notes: this code could be further optimized.
 * many things are saved in plain ascii instead of binary
 */

using namespace std;

#ifdef THREAD_SAFE_TIMING
static void atomic_double_add(std::atomic<double> &d1, double d2) {
      double current = d1.load();
        while (!d1.compare_exchange_weak(current, current + d2))
                ;
}
#endif
            
static atomic_double global_wtime_compactions (0), global_wtime_cdistribution (0), global_wtime_add_nodes (0), global_wtime_create_buckets (0), global_wtime_foreach_bucket (0), global_wtime_lambda (0), global_wtime_parallel (0), global_wtime_longest_lambda (0), global_wtime_best_sched(0);

static bool time_lambdas = true;
static std::mutex lambda_timing_mutex;
static size_t nb_threads_simulate=1; // this is somewhat a legacy parameter, i should get rid of (and replace by nb_threads)


namespace gatb { namespace core { namespace debruijn { namespace impl  {

    /* formerly lambda function inside bcalm but needed it in InsertIntoQueues also. no choice here  unless I wanted to typedef Model again*/
    #define minimizerMin(a,b) ((model.compareIntMinimizers(a,b)) ? a : b)
    #define minimizerMax(a,b) ((model.compareIntMinimizers(a,b)) ? b : a)

    /* class (formerly a simple lambda function) to process a kmer and decide which bucket(s) it should go to */
    /* needed to make it a class because i want it to remember its thread index */
    template <int SPAN>
    class InsertIntoQueues
    {
        typedef typename Kmer<SPAN>::Type  Type;
        typedef typename Kmer<SPAN>::Count Count;
        typedef typename Kmer<SPAN>::ModelCanonical ModelCanon;
        typedef typename Kmer<SPAN>::template ModelMinimizer <ModelCanon> Model;
        
        // new version, no longer using a queue-type object.
        typedef std::tuple<uint32_t, Type, uint32_t, uint32_t, uint32_t> tuple_t;
        typedef vector<tuple_t> flat_vector_queue_t;

        unsigned int p, k, abundance_threshold, nb_threads;
        Model &model, &modelK1;
        std::atomic<unsigned long>  &nb_left_min_diff_right_min, &nb_kmers_in_partition;
        Repartitor &repart;
        int _currentThreadIndex;
        std::vector<BankFasta*> &traveller_kmers_files;
        vector<std::mutex> &traveller_kmers_save_mutex;

        // saving traveller kmers in plain ASCII in files: a bit wasteful, but went to the easy solution
        void save_traveller_kmer (uint32_t minimizer, const string& seq, int abundance, uint32_t leftmin, uint32_t rightmin, int p) {
            Sequence s (Data::ASCII);
            s.getData().setRef ((char*)seq.c_str(), seq.size());
            s._comment = to_string(abundance); //abundance in comment
            traveller_kmers_save_mutex[p].lock();
            traveller_kmers_files[p]->insert(s);
            traveller_kmers_save_mutex[p].unlock();
        }

        public: 
        vector<flat_vector_queue_t> &flat_bucket_queues;

        /* function to add a kmer to a bucket */
        void add_to_bucket_queue(uint32_t minimizer, /*  string seq, */ Type &kmer, uint32_t abundance, uint32_t leftmin, uint32_t rightmin)
        {
            //bucket_queues.push_back(minimizer,std::make_tuple(TO_BUCKET_STR(seq),leftmin,rightmin,abundance));
            flat_bucket_queues[getThreadIndex()].push_back(std::make_tuple(minimizer, kmer, abundance, leftmin, rightmin));
        }

        /* boilerplate constructor */
        InsertIntoQueues(vector<flat_vector_queue_t> &flat_bucket_queues, 
                Model &model,
                Model &modelK1, 
                unsigned int p, unsigned int k, unsigned int nb_threads,
                int abundance_threshold,
                Repartitor &repart,
                std::atomic<unsigned long> &nb_left_min_diff_right_min,
                std::atomic<unsigned long> &nb_kmers_in_partition,
                std::vector<BankFasta*> &traveller_kmers_files,
                vector<std::mutex> &traveller_kmers_save_mutex
                ) : 
            p(p), k(k), abundance_threshold(abundance_threshold), nb_threads(nb_threads),
            model(model), modelK1(modelK1),
        nb_left_min_diff_right_min(nb_left_min_diff_right_min), nb_kmers_in_partition(nb_kmers_in_partition),
        repart(repart), _currentThreadIndex(-1), traveller_kmers_files(traveller_kmers_files),
        traveller_kmers_save_mutex(traveller_kmers_save_mutex),  flat_bucket_queues(flat_bucket_queues) {}

        /* does the actual work of processing a kmer, computing its minimizers, saving it to the right queue (basically the queue corresponding to its thread) */
        void operator()     (Count& item) {
            // if the abundance threshold is higher than the h5 abundance,
            // filter out this kmer (useful when you want to re-use same .h5 but with higher "-abundance" parameter)
            size_t abundance = item.abundance;
            if (abundance < (size_t)abundance_threshold)
                return;

            Type current = item.value; // current is a canonical kmer (i checked)
            uint32_t leftMin(modelK1.getMinimizerValue(current >> 2)); // that's because the lowest bit in the gatb kmer representation are the rightmost sequence nucleotides
            uint32_t rightMin(modelK1.getMinimizerValue(current));

            ++nb_kmers_in_partition;

            if (repart(leftMin) == p)
                add_to_bucket_queue(leftMin, current, abundance, leftMin, rightMin);

            if (leftMin != rightMin)
            {
                nb_left_min_diff_right_min ++;

                if (repart(rightMin) == p)
                    add_to_bucket_queue(rightMin, current, abundance, leftMin, rightMin);

                // handle "traveller kmers"
                uint32_t max_minimizer = minimizerMax(leftMin, rightMin);
                uint32_t min_minimizer = minimizerMin(leftMin, rightMin);
                if (repart(max_minimizer) != repart(min_minimizer))
                {
                    string seq = model.toString(current);
                    save_traveller_kmer(max_minimizer, seq, abundance, leftMin, rightMin, repart(max_minimizer));
                    //add_to_bucket_queue(max_minimizer, seq, leftMin, rightMin, repart(max_minimizer)); // no longer saved into the queue, but to a file instead

                    // sanity check
                    if (repart(max_minimizer) < repart(min_minimizer))
                    {                printf("unexpected problem: traveller kmer = %s, min_minimizer=%d max_minimizer=%d, repart(min_minimizer)=%d, repart(max_minimizer)=%d\n", seq.c_str(), min_minimizer, max_minimizer, repart(min_minimizer), repart(max_minimizer));                exit(1);            }
                }
            }

            // sanity check
            if (repart(leftMin) != p && repart(rightMin) != p)
            {                printf("unexpected problem: repart bucket\n");                exit(1);            }
        }

        /* neat trick taken from erwan's later work in gatb to find the thread id of a dispatched function */
        int getThreadIndex()
        {
            if (_currentThreadIndex < 0)
            {
                std::pair<IThread*,size_t> info;
                if (ThreadGroup::findThreadInfo (System::thread().getThreadSelf(), info) == true)
                {
                    _currentThreadIndex = info.second;
                }
                else
                {
                    throw Exception("Unable to find thread index during InsertIntoQueues");
                }
            }
            return _currentThreadIndex;
        }

    };

template<size_t SPAN>
void bcalm2(Storage *storage, 
        std::string prefix,
        int kmerSize, 
        int abundance_threshold, 
        int minSize, 
        int nb_threads, 
        int minimizer_type, 
        bool verbose
        )
{
    if (verbose)
        std::cout << "bcalm_algo params, prefix:" << prefix << " k:" << kmerSize << " a:" << abundance_threshold << " minsize:" << minSize << " threads:" << nb_threads << " mintype:" << minimizer_type << std::endl;
    if (nb_threads == 0) {std::cout << "oops bcalm called with 0 threads " << std::endl; exit(1);}
    if (minSize > kmerSize) {std::cout << "oops bcalm called with minSize(" << minSize << ")> kmerSize(" << kmerSize << ")" << std::endl; exit(1);}
    
    typedef typename Kmer<SPAN>::Type  Type;
    typedef typename Kmer<SPAN>::Count Count;
    typedef typename Kmer<SPAN>::ModelCanonical ModelCanon;
    typedef typename Kmer<SPAN>::template ModelMinimizer <ModelCanon> Model;

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);


    #ifdef BINSEQ
    initBinSeq(kmerSize);
    #endif
    
    /** We set BankBinary buffer. */
    BankBinary::setBufferSize (10000);

    auto start_t=chrono::system_clock::now();
    size_t maxBucket(0);

    /** We get the dsk and minimizers hash group in the storage object. */
    Group& dskGroup = storage->getGroup("dsk");
    Group& minimizersGroup = storage->getGroup("minimizers");

    typedef typename Kmer<SPAN>::Count Count;
    Partition<Count>& partition = dskGroup.getPartition<Count> ("solid");
    size_t nb_h5_partitions = partition.size();

    /* get actual number of partitions/passes _during_ DSK. oh so complicated .
     * this is needed because we need to group passes together. else 
     * the algo simply doesn't work. */
    Group& configGroup = storage->getGroup("configuration");
    stringstream ss; ss << configGroup.getProperty ("xml");
    Properties props; props.readXML (ss);
    size_t nb_passes = props.getInt("nb_passes");
    size_t nb_partitions = props.getInt("nb_partitions");

    if (verbose)
    {
        cout << "DSK used " << nb_passes << " passes and " << nb_partitions << " partitions" << std::endl;
    }

    if (nb_h5_partitions != nb_passes * nb_partitions)
    {
        cout << "Error: number of h5 partitions ("<< nb_h5_partitions << ") does not match number of DSK passes*partitions ("\
            << nb_passes<<"*"<<nb_partitions<<")" << endl;
        exit(1);
    }


    /** We retrieve the minimizers distribution from the solid kmers storage. */
    Repartitor repart;
    repart.load (minimizersGroup);

    u_int64_t rg = ((u_int64_t)1 << (2*minSize));

    /* Retrieve frequency of minimizers;
     * actually only used in minimizerMin and minimizerMax */
    uint32_t *freq_order = NULL;

    if (minimizer_type == 1)
    {
        freq_order = new uint32_t[rg];
        Storage::istream is (minimizersGroup, "minimFrequency");
        is.read ((char*)freq_order, sizeof(uint32_t) * rg);
    }

    Model model(kmerSize, minSize, typename Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);
    Model modelK1(kmerSize-1, minSize,  typename Kmer<SPAN>::ComparatorMinimizerFrequencyOrLex(), freq_order);

    std::vector<BankFasta*> out_to_glue(nb_threads); // each thread will write to its own glue file, to avoid locks
    
    // remove potential old glue files
    for (unsigned int i = 0; i < 10000 /* there cannot be more than 10000 threads, right? unsure if i'll pay for that asumption someday*/; i++)
    {
        if (System::file().doesExist(prefix + ".glue." + std::to_string(i)))
           System::file().remove (prefix + ".glue." + std::to_string(i)); 
    }

    unsigned long *nb_seqs_in_glue = new unsigned long[nb_threads];
    unsigned long *nb_pretips = new unsigned long[nb_threads];

    // another system could have been to send all sequences in a queue, and a thread responsible for writing to glue would dequeue (might be faster)
    for (unsigned int i = 0; i < (unsigned int)nb_threads; i++)
    {
        string glue_file = prefix + ".glue." + std::to_string(i);
        out_to_glue[i] = new BankFasta(glue_file);
        nb_seqs_in_glue[i] = 0;
        nb_pretips[i] = 0;
    }

    double weighted_best_theoretical_speedup_cumul = 0;
    double weighted_best_theoretical_speedup_sum_times = 0;
    double weighted_best_theoretical_speedup = 0;
    double weighted_actual_theoretical_speedup_cumul = 0;
    double weighted_actual_theoretical_speedup_sum_times = 0;
    double weighted_actual_theoretical_speedup = 0;

    auto start_buckets=chrono::system_clock::now();

    /* now our vocabulary is: a "DSK partition" == a "partition" == a "super-bucket" */
    /* buckets remain what they are in bcalm-original */
    /* a travelling kmer is one that goes to two buckets from different superbuckets */

    // I used to save traveller kmers into bucket_queues, but this would be a memory hog. Let's use files instead. Total volume will be small (a few gigs for human), but that's memory saved
    std::vector<BankFasta*> traveller_kmers_files(nb_partitions);
    vector<std::mutex> traveller_kmers_save_mutex(nb_partitions);
    std::string traveller_kmers_prefix = prefix + ".doubledKmers.";
    for (unsigned int i = 0; i < nb_partitions; i++)
        traveller_kmers_files[i] = new BankFasta(traveller_kmers_prefix + std::to_string(i));
   
    Dispatcher dispatcher (nb_threads); // setting up a multi-threaded dispatcher, so I guess we can say that things are getting pretty serious now

    // i want to do this but i'm not inside an Algorithm object:
    /*Iterator<int>* it_parts = Algorithm::createIterator<int>(
            new Range<int>::Iterator (0,nb_partitions-1), nb_partitions, "Iterating DSK partitions"
            );*/

    // copied from createIterator in Algorithm.hpp
    //  We create some listener to be notified every 1000 iterations and attach it to the iterator.
    IteratorListener* listener;
    if (verbose)
        listener = new ProgressTimer(nb_partitions, "Iterating DSK partitions");
    else
        listener = new IteratorListener ();

    auto it_parts = new tools::dp::impl::SubjectIterator<int> (
                new Range<int>::Iterator (0,nb_partitions-1), 
                nb_partitions/100);
    it_parts->addObserver (listener);
    LOCAL(it_parts);
    
    bcalm_logging = verbose;
    logging("prior to queues allocation");
 
    // new version, no longer using a queue-type object.
    typedef std::tuple<uint32_t, Type, uint32_t, uint32_t, uint32_t> tuple_t;
    typedef vector<tuple_t> flat_vector_queue_t;
    vector<flat_vector_queue_t> flat_bucket_queues(nb_threads);
       
    logging("Starting BCALM2");

    /*
     *
     * Iteration of partitions
     *
     *  main thread is going to read kmers from partitions and insert them into queues
     *
    */
    for (it_parts->first (); !it_parts->isDone(); it_parts->next()) /**FOREACH SUPERBUCKET (= partition) **/
    {
        uint32_t p = it_parts->item(); /* partition index */

        bool verbose_partition = verbose && ((p % ((nb_partitions+9)/10)) == 0); // only print verbose information 10 times at most

        size_t k = kmerSize;

        std::atomic<unsigned long> nb_left_min_diff_right_min;
        std::atomic<unsigned long> nb_kmers_in_partition;
        nb_kmers_in_partition = 0;
        nb_left_min_diff_right_min = 0;
        
        auto start_createbucket_t=get_wtime();
        
        InsertIntoQueues<SPAN> insertIntoQueues(flat_bucket_queues, model, modelK1, p, k, nb_threads, abundance_threshold, repart, nb_left_min_diff_right_min, nb_kmers_in_partition, traveller_kmers_files, traveller_kmers_save_mutex);

        /* MAIN FIRST LOOP: expand a superbucket by inserting kmers into queues. this creates buckets */
        // do it for all passes (because the union of passes correspond to a partition)
        for (size_t pass_index = 0 ; pass_index < nb_passes; pass_index ++)
        {
            /** We retrieve an iterator on the Count objects of the pth partition in pass pass_index */
            unsigned long interm_partition_index = p + pass_index * nb_partitions;
            Iterator<Count>* it_kmers = partition[interm_partition_index].iterator();
            LOCAL (it_kmers);

            if (pass_index == 0) // the first time, 
                for (int i = 0; i < nb_threads; i++) // resize approximately the bucket queues
                flat_bucket_queues[i].reserve(partition[interm_partition_index].getNbItems()/nb_threads);

            dispatcher.iterate (it_kmers, insertIntoQueues);
            /*for (it_kmers->first (); !it_kmers->isDone(); it_kmers->next()) // non-dispatcher version
                insertIntoQueues(it_kmers->item());*/
        }

        if (verbose_partition) 
            cout << endl << "Iterated " << nb_kmers_in_partition << " kmers, among them " << nb_left_min_diff_right_min << " were doubled" << endl;

        // also add traveller kmers that were saved to disk from a previous superbucket
        // but why don't we need to examine other partitions for potential traveller kmers?
        // no, because we iterate partitions in minimizer order.
        // but then you might say again something else:
        // "i thought bcalm1 needed to iterate partitions in minimizer order, but not bcalm2"
        // -> indeed, bcalm2 algorithm doesn't, but in the implementation i still choose to iterate in minimizer order.
        // because it seemed like a good idea at the time, when handling traveller kmers.
        // an alternative possibility would be to revert to minimizer-type 0 and repartition-type 0
        // advantages:
        // - this could enable loading multiple partitions at once (and more parallelization)
        // - faster kmer counting (16 mins vs 18 mins for cami medium, 1B distinct kmers)
        // but so far, I have not seen the need to load multiple partitions and the gain for dsk isnt big
        // disadvantages:
        // - would need to do a pass to write all traveller kmers to disk at first
        traveller_kmers_files[p]->flush();
        string traveller_kmers_file = traveller_kmers_prefix + std::to_string(p);
        std::atomic<unsigned long> nb_traveller_kmers_loaded;
        nb_traveller_kmers_loaded = 0;

        if (System::file().doesExist(traveller_kmers_file)) // for some partitions, there may be no traveller kmers
        {

            BankFasta traveller_kmers_bank (traveller_kmers_file);
            BankFasta::Iterator it (traveller_kmers_bank);
       
            class InsertTravellerKmer
            {
                int _currentThreadIndex;
                vector<flat_vector_queue_t> &flat_bucket_queues;
                Model &model, &modelK1;
                int k;
                std::atomic<unsigned long> &nb_traveller_kmers_loaded;

                public:
                InsertTravellerKmer(vector<flat_vector_queue_t> &flat_bucket_queues, Model& model, Model &modelK1, int k, std::atomic<unsigned long> &nb_traveller_kmers_loaded) 
                    : _currentThreadIndex(-1), flat_bucket_queues(flat_bucket_queues), model(model), modelK1(modelK1), k(k), nb_traveller_kmers_loaded(nb_traveller_kmers_loaded) {}

                int getThreadIndex()
                {
                    if (_currentThreadIndex < 0)
                    {
                        std::pair<IThread*,size_t> info;
                        if (ThreadGroup::findThreadInfo (System::thread().getThreadSelf(), info) == true)
                            _currentThreadIndex = info.second;
                        else
                            throw Exception("Unable to find thread index during InsertIntoQueues");
                    }
                    return _currentThreadIndex;
                }
                void operator () (const Sequence &sequence)
                {
                    string seq = sequence.toString();
                    string comment = sequence.getComment();
                    uint32_t abundance = atoi(comment.c_str());

                    // those could be saved in the BankFasta comment eventually
                    typename Model::Kmer current = model.codeSeed(seq.c_str(), Data::ASCII);
                    Type kmer = current.value();
                    uint32_t leftMin(modelK1.getMinimizerValue(kmer >> 2));
                    uint32_t rightMin(modelK1.getMinimizerValue(kmer));

                    uint32_t max_minimizer = minimizerMax(leftMin, rightMin);
                    //add_to_bucket_queue(max_minimizer, seq, abundance, leftMin, rightMin, p);
                    flat_bucket_queues[getThreadIndex()].push_back(std::make_tuple(max_minimizer, kmer, abundance, leftMin, rightMin));
                    nb_traveller_kmers_loaded++;
                }
            };
            InsertTravellerKmer insertTravellerKmer(flat_bucket_queues, model, modelK1, k, nb_traveller_kmers_loaded);

            dispatcher.iterate(it,insertTravellerKmer);

            if (verbose_partition) 
                std::cout << "Loaded " << nb_traveller_kmers_loaded << " doubled kmers for partition " << p << endl;
            traveller_kmers_bank.finalize();
            System::file().remove (traveller_kmers_file);
        }

        /* now that we have computed flat_bucket_queues' by each thread,
         * sort them by minimizer */

        //logging("begin sorting bucket queues");
        ThreadPool pool_sort(nb_threads);
        for (int thread = 0; thread < nb_threads; thread++)
        {
            auto sort_cmp = [] (tuple_t const &a, tuple_t const &b) -> bool { return get<0>(a) < get<0>(b); };

	    // todo check si  les minimiseurs sont pas deja quasiment triés dans un sens ou un autre, ca faciliterait le tri ici
            auto sort_bucket = [&sort_cmp, &flat_bucket_queues, thread] (int thread_id) 
            {std::sort(flat_bucket_queues[thread].begin(), flat_bucket_queues[thread].end(), sort_cmp);};

            if (nb_threads > 1)
                pool_sort.enqueue(sort_bucket);
            else
                sort_bucket(0);
        }
        pool_sort.join();
        //logging("end sorting bucket queues");

        /* remember which minimizer occurs in flat_bucket_queues' and its start position */
        set<uint32_t> set_minimizers;
        vector<uint64_t> nb_kmers_per_minimizer(rg);
        
        for (uint64_t i = 0; i < rg; i++)
            nb_kmers_per_minimizer[i] = 0;

        vector<vector<uint64_t>> start_minimizers(nb_threads);
        for (int thread = 0; thread < nb_threads; thread++)
        {
            // should be done in parallel possibly, if it takes time.
            set<uint32_t> set_minimizers_thread;
            start_minimizers[thread].resize(rg);
            uint64_t pos=0;
            //std:: cout << "iterating flat bucket queues  for thread " << thread << " elts: " << flat_bucket_queues[thread].size() << std::endl;
            for (auto v: flat_bucket_queues[thread])
            {
                uint32_t minimizer = get<0>(v);
                if (set_minimizers_thread.find(minimizer) == set_minimizers_thread.end())
                {
                    set_minimizers.insert(minimizer);
                    set_minimizers_thread.insert(minimizer);
                    start_minimizers[thread][minimizer] = pos;
                }
                nb_kmers_per_minimizer[minimizer]++;
                pos++;
            }
        }

        
        auto end_createbucket_t=get_wtime();
        atomic_double_add(global_wtime_create_buckets, diff_wtime(start_createbucket_t, end_createbucket_t));

        ThreadPool pool(nb_threads);

        std::vector<double> lambda_timings;
        auto start_foreach_bucket_t=get_wtime();

        /**FOREACH BUCKET **/
        for(auto actualMinimizer : set_minimizers)
        {
            auto lambdaCompact = [&nb_kmers_per_minimizer, actualMinimizer, &model,
                &maxBucket, &lambda_timings, &repart, &modelK1, &out_to_glue, &nb_seqs_in_glue, &nb_pretips, kmerSize, minSize,
                nb_threads, &start_minimizers, &flat_bucket_queues](int thread_id) {
                auto start_nodes_t=get_wtime();

                // (make sure to change other places labelled "// graph3" and "// graph4" as well)
                //graph4 g(kmerSize-1,actualMinimizer,minSize); // graph4
                uint number_elements(nb_kmers_per_minimizer[actualMinimizer]);
                #ifdef BINSEQ
                graph4 graphCompactor(kmerSize-1,actualMinimizer,minSize,number_elements);
                #else
                // cout<<"here"<<endl;
                //graph3 graphCompactor(kmerSize-1,actualMinimizer,minSize,number_elements);
                graph3<SPAN> graphCompactor(kmerSize-1,actualMinimizer,minSize,number_elements); // graph3<span> switch 
                graphCompactor.pre_tip_cleaning = false; // this is the actual trigger for bcalm pre-tip simplifications. 
                                                        // i'm leaving it off for now because the gains do not seem that big
                #endif

                /* add nodes to graph */
                //while (bucket_queues.pop_immediately(actualMinimizer,bucket_elt))

                /* go through all the flat_bucket_queues's that were constructed by each thread,
                     * and iterate a certain minimizer. i dont even need a priority queue! */

                // used to be in a lambda outside of that lambda, there was a bug, decided to put it here but didnt even solve the bug, hmm. i should have been more explicit whether the bug still happens or not, i dunno now.
                for (int thread = 0; thread < nb_threads; thread++)
                {
                    uint64_t pos = start_minimizers[thread][actualMinimizer];
                    unsigned int size = flat_bucket_queues[thread].size();
                    if (pos == size) continue;
                    while (actualMinimizer == get<0>(flat_bucket_queues[thread][pos]))
                    {
                        auto tupl = flat_bucket_queues[thread][pos]; // the tuple format in flat_bucket_queues is: (minimizer, seq, abundance, leftmin, rightmin)
                        std::tuple<BUCKET_STR_TYPE,uint,uint,uint> bucket_elt; // graph3<span> switch 
                        // g.addleftmin(std::get<1>(bucket_elt));
                        // g.addrightmin(std::get<2>(bucket_elt));
                        // g.addvertex(FROM_BUCKET_STR(std::get<0>(bucket_elt)));
                        string seq = model.toString(get<1>(tupl));
                        uint32_t a = get<2>(tupl), b = get<3>(tupl), c = get<4>(tupl);
                        bucket_elt = make_tuple(seq,b,c,a);
                        //std::cout << " (debug) adding to graph: " << std::get<0>(bucket_elt) << std::endl;
                        graphCompactor.addtuple(bucket_elt); // addtuple wants that tuple: (seq, leftmin, rightmin, abundance)

                        pos++;
                        if (pos == size) break;
                    }
                }


                // cout<<"endaddtuple"<<endl;
                auto end_nodes_t=get_wtime();
                atomic_double_add(global_wtime_add_nodes, diff_wtime(start_nodes_t, end_nodes_t));

                /* compact graph*/
                auto start_dbg=get_wtime();
                graphCompactor.debruijn();

                auto end_dbg=get_wtime();
                atomic_double_add(global_wtime_compactions, diff_wtime(start_dbg, end_dbg));

                /* distribute nodes (to other buckets, or output, or glue) */
                auto start_cdistribution_t=get_wtime();
                string seq;
                for(uint32_t i(0);i<number_elements;++i){
                    if(graphCompactor.output(i)){
                        #ifdef BINSEQ
    					seq=graphCompactor.unitigs[i].str(); // graph4
                        #else
                        seq=graphCompactor.unitigs[i]; // graph3
                        //std::vector<unsigned int> abundances ; // graph3 
                        std::vector<unsigned int>& abundances = graphCompactor.unitigs_abundances[i]; // graph3 // graph3<span> switch 
                        #endif
                        //std::cout << " (debug) got from compacted graph: " << seq << std::endl;

                        typename Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, kmerSize - 1).c_str(), Data::ASCII);
                        uint leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                        typename Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - kmerSize + 1, kmerSize - 1).c_str(), Data::ASCII);
                        uint rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
                        bool lmark = actualMinimizer != leftMin;
                        bool rmark = actualMinimizer != rightMin;

                        Sequence s (Data::ASCII);
                        s.getData().setRef ((char*)seq.c_str(), seq.size());
                        s._comment = string(lmark?"1":"0")+string(rmark?"1":"0"); //We set the sequence comment.
                        s._comment += " ";
                        for (auto abundance : abundances)
                            s._comment += to_string(abundance) + " ";
                        out_to_glue[thread_id]->insert(s); 
                        nb_seqs_in_glue[thread_id]++;
                    }
                }
                nb_pretips[thread_id] += graphCompactor.nb_pretips;
                graphCompactor.nb_pretips = 0;
                graphCompactor.clear(); // frees memory allocated during graph3 constructor (sort of a destructor, if you will)
                auto end_cdistribution_t=get_wtime();
                atomic_double_add(global_wtime_cdistribution, diff_wtime(start_cdistribution_t, end_cdistribution_t));

                if(number_elements>maxBucket){maxBucket=number_elements;}

                if (time_lambdas)
                {
                    auto time_lambda = diff_wtime(start_nodes_t, end_cdistribution_t);
                    atomic_double_add(global_wtime_lambda, time_lambda);
                    lambda_timing_mutex.lock();
                    lambda_timings.push_back(time_lambda);
                    lambda_timing_mutex.unlock();
                }

            }; // end lambda function

            if (nb_threads > 1)
                pool.enqueue(lambdaCompact);
            else
                lambdaCompact(0);

        } // end for each bucket

        pool.join();
        //logging("done compactions");
            
        // flush glues, clear flat_bucket_queues
        for (int thread_id = 0; thread_id < nb_threads; thread_id++)
        {
            flat_bucket_queues[thread_id].clear();
            out_to_glue[thread_id]->flush (); 
        }

        if (partition[p].getNbItems() == 0)
            continue; // no stats to print here

        /* compute and print timings */
        {
            auto end_foreach_bucket_t=get_wtime();
            auto wallclock_sb = diff_wtime(start_foreach_bucket_t, end_foreach_bucket_t);
            atomic_double_add(global_wtime_foreach_bucket, wallclock_sb);
            atomic_double_add(global_wtime_parallel, wallclock_sb);

            if (time_lambdas && lambda_timings.size() > 0)
            {
                std::sort(lambda_timings.begin(), lambda_timings.end());
                std::reverse(lambda_timings.begin(), lambda_timings.end());
                /* compute a theoretical, i think optimal, scheduling of lambda's using the current number of threads
                */
                double tot_time_best_sched_lambda = 0; // start with the longest lambda
                int t = 0;
                for (auto & lambda_time: lambda_timings) {
                    if ((t++) % nb_threads_simulate == 0)
                        tot_time_best_sched_lambda += lambda_time;
                }

                double longest_lambda = lambda_timings.front();

                if (verbose_partition)
                {
                    cout <<"\nIn this superbucket (containing " << set_minimizers.size() << " active minimizers)," <<endl;
                    cout <<"                  sum of time spent in lambda's: "<< global_wtime_lambda / 1000000 <<" msecs" <<endl;
                    cout <<"                                 longest lambda: "<< longest_lambda / 1000000 <<" msecs" <<endl;
                    cout <<"         tot time of best scheduling of lambdas: "<< tot_time_best_sched_lambda / 1000000 <<" msecs" <<endl;
                }

                double best_theoretical_speedup =  global_wtime_lambda  / longest_lambda;
                double actual_theoretical_speedup =  global_wtime_lambda  / tot_time_best_sched_lambda;

                if (verbose_partition)
                {
                    cout <<"                       best theoretical speedup: "<<  best_theoretical_speedup << "x" <<endl;
                    if (nb_threads_simulate > 1)
                        cout <<"     best theoretical speedup with "<< nb_threads_simulate << " thread(s): "<<  actual_theoretical_speedup << "x" <<endl;
                }

                weighted_best_theoretical_speedup_cumul += best_theoretical_speedup * wallclock_sb;
                weighted_best_theoretical_speedup_sum_times                        += wallclock_sb;
                weighted_best_theoretical_speedup = weighted_best_theoretical_speedup_cumul / weighted_best_theoretical_speedup_sum_times ;

                weighted_actual_theoretical_speedup_cumul += actual_theoretical_speedup * wallclock_sb;
                weighted_actual_theoretical_speedup_sum_times                        += wallclock_sb;
                weighted_actual_theoretical_speedup = weighted_actual_theoretical_speedup_cumul / weighted_actual_theoretical_speedup_sum_times ;
                atomic_double_add(global_wtime_longest_lambda, longest_lambda);
                atomic_double_add(global_wtime_best_sched, tot_time_best_sched_lambda);
                global_wtime_lambda = 0;
            }
        }

        if (verbose_partition)
            logging("Done with partition " + std::to_string(p));
    } // end iteration superbuckets
    
    /*
     *
     * Finishing up
     *
     */

    // create list of non empty glues
    std::ofstream list_of_glues(prefix + ".glue"); 
    for (unsigned int i = 0; i < (unsigned int)nb_threads; i++)
    {

        char* prefix_copy = strdup (prefix.c_str());
        std::string prefix_base (basename(prefix_copy)); // posix basename() may alter the string
        free (prefix_copy);

        string glue_file = prefix_base + ".glue." + std::to_string(i);
        if (nb_seqs_in_glue[i])
        {
            list_of_glues << glue_file << endl;
        }
    }
    list_of_glues.close();

    // gather some stats
    uint64_t nbSeqsInGlue = 0;
    uint64_t nbPretips = 0;
    for (int thread_id = 0; thread_id < nb_threads; thread_id++)
    {
       nbSeqsInGlue += nb_seqs_in_glue[thread_id];
       nbPretips += nb_pretips[thread_id];
    }

    auto end_t=chrono::system_clock::now();
    float wtime = chrono::duration_cast<chrono::nanoseconds>(end_t - start_buckets).count() / unit;

    Group& bcalmGroup = storage->getGroup("bcalm");
    bcalmGroup.setProperty ("nb_pretips_removed",     Stringify::format("%ld", nb_pretips));
    bcalmGroup.setProperty ("nb_sequences_in_glue",     Stringify::format("%ld", nbSeqsInGlue));
    bcalmGroup.setProperty ("wtime_compactions",     Stringify::format("%f", wtime));

    /* printing some timing stats */
    if (verbose) 
    {
        cout <<"Number of sequences in glue: "<< nbSeqsInGlue << std::endl;
        cout <<"Number of pre-tips removed : "<< nbPretips << std::endl;
        cout<<"Buckets compaction and gluing           : "<< wtime <<" secs"<<endl;
        cout<<"Within that, \n";
        cout <<"                                 creating buckets from superbuckets: "<< global_wtime_create_buckets / unit <<" secs"<<endl;
        cout <<"                      bucket compaction (wall-clock during threads): "<< global_wtime_foreach_bucket / unit <<" secs" <<endl;
        cout <<"\n                within all bucket compaction threads,\n";
        cout <<"                       adding nodes to subgraphs: "<< global_wtime_add_nodes / unit <<" secs" <<endl;
        cout <<"         subgraphs constructions and compactions: "<< global_wtime_compactions / unit <<" secs"<<endl;
        cout <<"                  compacted nodes redistribution: "<< global_wtime_cdistribution / unit <<" secs"<<endl;
    }
    double sum =  global_wtime_cdistribution + global_wtime_compactions + global_wtime_add_nodes + global_wtime_create_buckets;
    if (verbose) 
        cout<<"Sum of CPU times for bucket compactions: "<< sum / unit <<" secs"<<endl;
    if (nb_threads == 1)
    {
        if (verbose)
        cout<<"Discrepancy between sum of fine-grained timings and total wallclock of buckets compactions step: "<< (chrono::duration_cast<chrono::nanoseconds>(end_t-start_buckets).count() - sum ) / unit <<" secs"<<endl;
    }

    if (verbose)
        cout<<"BCALM total wallclock (excl kmer counting): "<<chrono::duration_cast<chrono::nanoseconds>(end_t-start_t).count() / unit <<" secs"<<endl;

    if (verbose) 
        cout<<"Maximum number of kmers in a subgraph: "<<maxBucket<<endl;

    if (time_lambdas && verbose)
    {
        cout<<"Performance of compaction step:\n"<<endl;
        cout<<"                 Wallclock time spent in parallel section : "<< global_wtime_parallel / unit << " secs"<<endl;
        cout<<"             Best theoretical speedup in parallel section : "<< weighted_best_theoretical_speedup << "x" <<endl;
        cout<<"Best theor. speedup in parallel section using " << nb_threads_simulate << " threads : "<< weighted_actual_theoretical_speedup << "x" <<endl;
        cout<<"             Sum of longest bucket compaction for each sb : "<< global_wtime_longest_lambda / unit << " secs"<<endl;
        cout<<"                       Sum of best scheduling for each sb : "<< global_wtime_best_sched / unit << " secs"<<endl;
    }

    // cleanup everything that was new'd
    if (minimizer_type == 1)
        delete[] freq_order;
    delete[] nb_seqs_in_glue;
    delete[] nb_pretips;
    for (unsigned int i = 0; i < (unsigned int)nb_threads; i++)
        delete out_to_glue[i];
    for (unsigned int i = 0; i < nb_partitions; i++)
        delete traveller_kmers_files[i];

    logging("Done with all compactions");

    //delete storage; exit(0); // to stop after bcalm, before bglue
}


}}}}
