#include "bcalm_algo.hpp"
#include <libgen.h> // for basename()

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
static std::mutex lambda_timing_mutex, active_minimizers_mutex;
static size_t nb_threads_simulate=1; // this is somewhat a legacy parameter, i should get rid of (and replace by nb_threads)


static unsigned long memory_usage(string message="", bool verbose=true)
{
    // using Progress.cpp of gatb-core
    u_int64_t mem = System::info().getMemorySelfUsed() / 1024;
    u_int64_t memMaxProcess = System::info().getMemorySelfMaxUsed() / 1024;
    char tmp[128];
    snprintf (tmp, sizeof(tmp), "  --  memory [current, maximum (maxRSS)]: [%4lu, %4lu] MB ",
            mem, memMaxProcess);
    if (verbose)
    {        
        std::cout << message << " " << tmp << std::endl;
    }
    return mem;
}


namespace gatb { namespace core { namespace debruijn { namespace impl  {

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

    auto minimizerMin = [&repart, &model] (uint32_t a, uint32_t b)
    {
        return (model.compareIntMinimizers(a,b)) ? a : b;
    };

    auto minimizerMax = [&repart, &model] (uint32_t a, uint32_t b)
    {
        return (model.compareIntMinimizers(a,b)) ? b : a;
    };

    std::vector<BankFasta*> out_to_glue(nb_threads); // each thread will write to its own glue file, to avoid locks
    
    // remove potential old glue files
    for (unsigned int i = 0; i < 10000 /* there cannot be more than 10000 threads, right? unsure if i'll pay for that asumption someday*/; i++)
        if (System::file().doesExist(prefix + ".glue." + std::to_string(i)))
        {
           System::file().remove (prefix + ".glue." + std::to_string(i)); 
        }

    unsigned long *nb_seqs_in_glue = new unsigned long[nb_threads];

    // another system could have been to send all sequences in a queue, and a thread responsible for writing to glue would dequeue (might be faster)
    for (unsigned int i = 0; i < (unsigned int)nb_threads; i++)
    {
        string glue_file = prefix + ".glue." + std::to_string(i);
        out_to_glue[i] = new BankFasta(glue_file);
        nb_seqs_in_glue[i] = 0;
    }

    double weighted_best_theoretical_speedup_cumul = 0;
    double weighted_best_theoretical_speedup_sum_times = 0;
    double weighted_best_theoretical_speedup = 0;
    double weighted_actual_theoretical_speedup_cumul = 0;
    double weighted_actual_theoretical_speedup_sum_times = 0;
    double weighted_actual_theoretical_speedup = 0;

    auto start_buckets=chrono::system_clock::now();

    std::vector<std::set<uint32_t>> active_minimizers;
    active_minimizers.resize(nb_partitions);

     /* now our vocabulary is: a "DSK partition" == a "partition" == a "super-bucket" */
    /* buckets remain what they are in bcalm-original */
    /* a travelling kmer is one that goes to two buckets from different superbuckets */

    // I used to save traveller kmers into bucket_queues, but this would be a memory hog. Let's use files instead. Total volume will be small (a few gigs for human), but that's memory saved
    std::vector<BankFasta*> traveller_kmers_files(nb_partitions);
    std::string traveller_kmers_prefix = prefix + ".doubledKmers.";
    std::mutex *traveller_kmers_save_mutex = new std::mutex[nb_partitions];
    for (unsigned int i = 0; i < nb_partitions; i++)
        traveller_kmers_files[i] = new BankFasta(traveller_kmers_prefix + std::to_string(i));

    auto save_traveller_kmer = [&traveller_kmers_files, &traveller_kmers_save_mutex]
        (uint32_t minimizer, string seq, int abundance, uint32_t leftmin, uint32_t rightmin, int p) {
            // saving traveller kmers in plain ASCII in files: a bit wasteful, but went to the easy solution
            Sequence s (Data::ASCII);
            s.getData().setRef ((char*)seq.c_str(), seq.size());
            s._comment = to_string(abundance); //abundance in comment
            traveller_kmers_save_mutex[p].lock();
            traveller_kmers_files[p]->insert(s);
            traveller_kmers_files[p]->flush();
            traveller_kmers_save_mutex[p].unlock();

        };
    
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

        size_t k = kmerSize;

        // create many queues in place of Buckets
        // (this code used to be outside the partition loop, but I think it's a good idea to reinit the queues after each superbucket(=partition) to avoid queues leaking memory

        // this implementation is supposedly efficient, but:
        // - as fast as the lockbasedqueue below
        // - uses much more memory
        //moodycamel::ConcurrentQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> > bucket_queues[rg];

        // another queue system, very simple, with locks
        // it's fine but uses a linked list, so more memory than I'd like
        //LockBasedQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> > bucket_queues[rg];

        // still uses more memory than i'd like
        // LockStdQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> > bucket_queues[rg];

        //LockStdQueue<std::tuple<BUCKET_STR_TYPE,uint32_t, uint32_t> > *bucket_queues=new LockStdQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t > > [rg];
        LockStdQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t, uint32_t> > *bucket_queues=new LockStdQueue<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t, uint32_t> > [rg]; // graph3<span> switch 

        //LockStdVector<std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> > bucket_queues[rg]; // very inefficient


        /* lambda function to add a kmer to a bucket */
        auto add_to_bucket_queue = [&active_minimizers, &bucket_queues](uint32_t minimizer, string seq, int abundance, uint32_t leftmin, uint32_t rightmin, int p)
        {
            //std::cout << "adding elt to bucket: " << seq << " "<< minimizer<<std::endl;
            //bucket_queues[minimizer].enqueue(std::make_tuple(TO_BUCKET_STR(seq),leftmin,rightmin));
            bucket_queues[minimizer].enqueue(std::make_tuple(TO_BUCKET_STR(seq),leftmin,rightmin,abundance)); // graph3<span> switch

            if (active_minimizers[p].find(minimizer) == active_minimizers[p].end())
            {
                active_minimizers_mutex.lock();
                active_minimizers[p].insert(minimizer);
                active_minimizers_mutex.unlock();
            }
        };

        std::atomic<unsigned long> nb_left_min_diff_right_min;
        std::atomic<unsigned long> nb_kmers_in_partition;
        nb_kmers_in_partition = 0;
        nb_left_min_diff_right_min = 0;
        std::atomic<uint32_t> kmerInGraph;
        kmerInGraph = 0;

        /* lambda function to process a kmer and decide which bucket(s) it should go to */
        auto insertIntoQueues = [p, &minimizerMax, &minimizerMin, &add_to_bucket_queue,
                    &bucket_queues, &modelK1, &k, &repart, &nb_left_min_diff_right_min,
                    &kmerInGraph, &model, &save_traveller_kmer, abundance_threshold,
                    &nb_kmers_in_partition]
                (Count& item) {
            
            // if the abundance threshold is higher than the h5 abundance,
            // filter out this kmer (useful when you want to re-use same .h5 but with higher "-abundance" parameter)
            size_t abundance = item.abundance;
            if (abundance < (size_t)abundance_threshold)
                return;

            Type current = item.value;

            string seq = model.toString(current);
            typename Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
            uint32_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
            typename Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
            uint32_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
            // string seq;
            // uint32_t leftMin(0);
            // uint32_t rightMin(0);

            ++kmerInGraph;
            ++nb_kmers_in_partition;
            
            if (repart(leftMin) == p)
                add_to_bucket_queue(leftMin, seq, abundance, leftMin, rightMin, p);

            if (leftMin != rightMin)
            {
                nb_left_min_diff_right_min ++;

                if (repart(rightMin) == p)
                    add_to_bucket_queue(rightMin, seq, abundance, leftMin, rightMin, p);

                // handle traveller kmers
                uint32_t max_minimizer = minimizerMax(leftMin, rightMin);
                uint32_t min_minimizer = minimizerMin(leftMin, rightMin);
                if (repart(max_minimizer) != repart(min_minimizer))
                {
                    /* I call that a "traveller kmer" */
                    save_traveller_kmer(max_minimizer, seq, abundance, leftMin, rightMin, repart(max_minimizer));
                    //add_to_bucket_queue(max_minimizer, seq, leftMin, rightMin, repart(max_minimizer)); // no longer saved into the queue, but to a file instead

                    // sanity check
                    if (repart(max_minimizer) < repart(min_minimizer))
                    {                printf("wtf? traveller kmer = %s, min_minimizer=%d max_minimizer=%d, repart(min_minimizer)=%d, repart(max_minimizer)=%d\n", seq.c_str(), min_minimizer, max_minimizer, repart(min_minimizer), repart(max_minimizer));                exit(1);            }
                }
            }

            // sanity check
            if (repart(leftMin) != p && repart(rightMin) != p)
            {                printf("wtf? repart bucket\n");                exit(1);            }

        };
        

        auto start_createbucket_t=get_wtime();

        /* MAIN FIRST LOOP: expand a superbucket by inserting kmers into queues. this creates buckets */
        // do it for all passes (because the union of passes correspond to a partition)
        for (size_t pass_index = 0 ; pass_index < nb_passes; pass_index ++)
        {
            /** We retrieve an iterator on the Count objects of the pth partition in pass pass_index */
            unsigned long interm_partition_index = p + pass_index * nb_partitions;
            Iterator<Count>* it_kmers = partition[interm_partition_index].iterator();
            LOCAL (it_kmers);
            dispatcher.iterate (it_kmers, insertIntoQueues);
        }

        if (verbose) 
            cout << endl << "Iterated " << nb_kmers_in_partition << " kmers, among them " << nb_left_min_diff_right_min << " were doubled" << endl;

        // also add traveller kmers that were saved to disk from a previous superbucket
        // but why don't we need to examine other partitions for potential traveller kmers?
        // no, because we iterate partitions in minimizer order.
        // but then you might again something else:
        // "i thought bcalm1 needed to iterate partitions in minimizer order, but not bcalm2"
        // -> indeed, bcalm2 algorithm doesn't, but in the implementation i still choose to iterate in minimizer order.
        // because it seemed like a good idea at the time, when handling traveller kmers.
        // looking back, it might be a good idea to not do that anymore.
        // this could enable loading multiple partitions at once (and more parallelization)
        string traveller_kmers_file = traveller_kmers_prefix + std::to_string(p);
        unsigned long nb_traveller_kmers_loaded = 0;
        if (System::file().doesExist(traveller_kmers_file)) // for some partitions, there may be no traveller kmers
        {
            BankFasta traveller_kmers_bank (traveller_kmers_file);
            BankFasta::Iterator it (traveller_kmers_bank);
            for (it.first(); !it.isDone(); it.next())
            {
                string seq = it->toString();
                string comment = it->getComment();
                int abundance = atoi(comment.c_str());

                // those could be saved in the BankFasta comment eventually
                typename Model::Kmer kmmerBegin = modelK1.codeSeed(seq.substr(0, k - 1).c_str(), Data::ASCII);
                uint32_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
                typename Model::Kmer kmmerEnd = modelK1.codeSeed(seq.substr(seq.size() - k + 1, k - 1).c_str(), Data::ASCII);
                uint32_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));

                uint32_t max_minimizer = minimizerMax(leftMin, rightMin);
                add_to_bucket_queue(max_minimizer, seq, abundance, leftMin, rightMin, p);
                nb_traveller_kmers_loaded++;
            }
            if (verbose) 
                std::cout << "Loaded " << nb_traveller_kmers_loaded << " doubled kmers for partition " << p << endl;
            traveller_kmers_bank.finalize();
            System::file().remove (traveller_kmers_file);
        }
        

        auto end_createbucket_t=get_wtime();
        atomic_double_add(global_wtime_create_buckets, diff_wtime(start_createbucket_t, end_createbucket_t));

        ThreadPool pool(nb_threads);

        std::vector<double> lambda_timings;
        auto start_foreach_bucket_t=get_wtime();

        /**FOREACH BUCKET **/
        for(auto actualMinimizer : active_minimizers[p])
        {
            auto lambdaCompact = [&bucket_queues, actualMinimizer,
                &maxBucket, &lambda_timings, &repart, &modelK1, &out_to_glue, &nb_seqs_in_glue, kmerSize, minSize](int thread_id) {
                auto start_nodes_t=get_wtime();

                bool debug = false;
                
                // (make sure to change other places labelled "// graph3" and "// graph4" as well)
                //graph4 g(kmerSize-1,actualMinimizer,minSize); // graph4
                uint number_elements(bucket_queues[actualMinimizer].size_approx());
                #ifdef BINSEQ
                graph4 graphCompactor(kmerSize-1,actualMinimizer,minSize,number_elements);
                #else
                // cout<<"here"<<endl;
                //graph3 graphCompactor(kmerSize-1,actualMinimizer,minSize,number_elements);
                graph3<SPAN> graphCompactor(kmerSize-1,actualMinimizer,minSize,number_elements); // graph3<span> switch 
                #endif

                /* add nodes to graph */
                //std::tuple<BUCKET_STR_TYPE,uint,uint> bucket_elt;
                std::tuple<BUCKET_STR_TYPE,uint,uint,uint> bucket_elt; // graph3<span> switch 

                while (bucket_queues[actualMinimizer].try_dequeue(bucket_elt))
                {
                // for(uint i(0);i<number_elements;++i)
                // {
                //     bucket_queues[actualMinimizer].try_dequeue(bucket_elt);
                    // g.addleftmin(std::get<1>(bucket_elt));
                    // g.addrightmin(std::get<2>(bucket_elt));
                    // g.addvertex(FROM_BUCKET_STR(std::get<0>(bucket_elt)));
                    if (debug) 
                        std::cout << " (debug) adding to graph: " << std::get<0>(bucket_elt) << std::endl;
                    graphCompactor.addtuple(bucket_elt);
                   
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
                        if (debug)
                        std::cout << " (debug) got from compacted graph: " << seq << std::endl;

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

        // flush glues
        for (unsigned int thread_id = 0; thread_id < (unsigned int)nb_threads; thread_id++)
        {
            out_to_glue[thread_id]->flush (); 
        }


        if (partition[p].getNbItems() == 0)
            continue; // no stats to print here

        // check if buckets are indeed empty
        for (unsigned int minimizer = 0; minimizer < rg; minimizer++)
        {
            if  (bucket_queues[minimizer].size_approx() != 0)
            {
                printf("WARNING! bucket %d still has non-processed %d elements\n", minimizer, bucket_queues[minimizer].size_approx() );
                //std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t> bucket_elt; 
                std::tuple<BUCKET_STR_TYPE,uint32_t,uint32_t,uint32_t> bucket_elt; // graph3<span> switch 
                while (bucket_queues[minimizer].try_dequeue(bucket_elt))
                {
                    printf("    %s leftmin %d rightmin %d abundance %d repartleft %d repartright %d repartmin %d\n", FROM_BUCKET_STR(std::get<0>(bucket_elt)).c_str(), std::get<1>(bucket_elt), std::get<2>(bucket_elt), std::get<3>(bucket_elt), repart(std::get<1>(bucket_elt)), repart(std::get<2>(bucket_elt)), repart(minimizer)); // graph3<span> switch 
                }

            }
        }


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

                if (verbose) 
                {
                    cout <<"\nIn this superbucket (containing " << active_minimizers.size() << " active minimizers)," <<endl;
                    cout <<"                  sum of time spent in lambda's: "<< global_wtime_lambda / 1000000 <<" msecs" <<endl;
                    cout <<"                                 longest lambda: "<< longest_lambda / 1000000 <<" msecs" <<endl;
                    cout <<"         tot time of best scheduling of lambdas: "<< tot_time_best_sched_lambda / 1000000 <<" msecs" <<endl;
                }

                double best_theoretical_speedup =  global_wtime_lambda  / longest_lambda;
                double actual_theoretical_speedup =  global_wtime_lambda  / tot_time_best_sched_lambda;

                if (verbose) 
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

        delete [] bucket_queues;
        memory_usage("Done with partition " + std::to_string(p), verbose);
    } // end iteration superbuckets


        // FIXME there may be a memory leak here, test it (saw it on spruce)

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

    /* printing some timing stats */
    auto end_t=chrono::system_clock::now();
    if (verbose) 
    {
        cout<<"Buckets compaction and gluing           : "<<chrono::duration_cast<chrono::nanoseconds>(end_t - start_buckets).count() / unit<<" secs"<<endl;
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
    for (unsigned int i = 0; i < (unsigned int)nb_threads; i++)
        delete out_to_glue[i];
    delete[] traveller_kmers_save_mutex;
    for (unsigned int i = 0; i < nb_partitions; i++)
        delete traveller_kmers_files[i];
 
    memory_usage("Done with all compactions", verbose);

    //delete storage; exit(0); // to stop after bcalm, before bglue
}

}}}}
