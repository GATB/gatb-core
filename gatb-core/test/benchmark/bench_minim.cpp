
#include <chrono>
#define get_wtime() chrono::system_clock::now()
#define diff_wtime(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count()


#include <gatb/system/impl/System.hpp>

#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/Traversal.hpp>

#include <gatb/bank/impl/BankStrings.hpp>
#include <gatb/bank/impl/BankSplitter.hpp>
#include <gatb/bank/impl/BankRandom.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <gatb/kmer/impl/Model.hpp>


#include <iostream>
#include <memory>

using namespace std;

using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;



#define span KMER_SPAN(0) // fixme this need to be modified for larger k's



/* inspired by debruijn_test3 from unit tests*/


void debruijn_minim_bench (const Graph& graph, size_t kmerSize)
{
    int miniSize = 8;
    int NB_REPETITIONS = 2000000;

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(3);

    Graph::Iterator<Node> nodes = graph.iterator<Node> ();
    nodes.first ();

    /** We get the first node. */
    Node node = nodes.item();

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::ModelCanonical  ModelCanonical;
    typedef typename Kmer<span>::ModelDirect     ModelDirect;
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   ModelMini;
    typedef typename ModelMini::Kmer                        KmerType;

    ModelMini  modelMini (kmerSize, miniSize);
    
    if (span == KMER_SPAN(1)) 
        cout << "this should be printed right before the segfault" << endl;

    /** We get the value of the current minimizer. */
    Type kmer = node.kmer.get<Type>();

    if (span == KMER_SPAN(1)) 
        cout << "this should be printed right after the segfault" << endl;

    auto start_t=chrono::system_clock::now();
    for (unsigned int i = 0 ; i < NB_REPETITIONS ; i++)
        modelMini.getMinimizerValue(kmer);
    auto end_t=chrono::system_clock::now();
    cout << "start of function 4" << endl;
			
    cout << NB_REPETITIONS << " minimizers of length " << miniSize << " on a " << kmerSize << "-mer : " << diff_wtime(start_t, end_t) / unit << " seconds" << endl;
	
    cout << "---- now on all nodes of the graph -----\n";

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
    {
        /* do nothing */
    }
    end_t=chrono::system_clock::now();
    auto baseline_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline iteration time (" << nodes.size() << " nodes) : " << baseline_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    int times = max((int)(NB_REPETITIONS / nodes.size()), 1);
    for (int i=0; i < times; i++)
        for (nodes.first(); !nodes.isDone(); nodes.next())
            modelMini.getMinimizerValue(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();
    cout << nodes.size() << " minimizers of length " << miniSize << " on all nodes (" << kmerSize << "-mers), " << to_string(times) << " times, with existing code : " << (diff_wtime(start_t, end_t) / unit) - baseline_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (int i=0; i < times; i++)
        for (nodes.first(); !nodes.isDone(); nodes.next())
            modelMini.getMinimizerValueNew(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();
    cout << nodes.size() << " minimizers of length " << miniSize << " on all nodes (" << kmerSize << "-mers), " << to_string(times) << " times, with new method    : " << (diff_wtime(start_t, end_t) / unit) - baseline_time << " seconds" << endl;

    
    cout << "checking agreement... ";
    for (nodes.first(); !nodes.isDone(); nodes.next())
    {
        if (modelMini.getMinimizerValue(nodes.item().kmer.get<Type>()) != modelMini.getMinimizerValueNew(nodes.item().kmer.get<Type>()))
        {
            cout << "FAIL! problem with kmer " << graph.toString(node) << " : " << modelMini.getMinimizerValue(nodes.item().kmer.get<Type>()) << " vs " << modelMini.getMinimizerValueNew(nodes.item().kmer.get<Type>()) << endl;
            exit(1);
        }
    }
    cout << "all good." << endl;
}


void debruijn_minim ()
{
    const char* sequences [] =
    {
//        "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAAT",
//        "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA",
        "CGCTATTCATCATTGTTTATCAATGAGCTAAAAGGAAACTATAAATAACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
    };

    //        size_t kmerSizes[] = { 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};
    size_t kmerSizes[] = {21, 31};

    for (size_t i=0; i<ARRAY_SIZE(sequences); i++)
    {
        for (size_t j=0; j<ARRAY_SIZE(kmerSizes); j++)
        {
            /** We create the graph. */
            Graph graph = Graph::create (
                    new BankStrings (sequences[i], 0),
                    "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d -mphf emphf", kmerSizes[j], 500 
                    );

            debruijn_minim_bench (graph, kmerSizes[j]);

            /** We remove the graph. */
            graph.remove ();
        }
    }
}

int main (int argc, char* argv[])
{
    try
    {
        // if no arg provided, just run the basic test on a single node
        if (argc == 1)
            debruijn_minim();
        else
            // else, use a real file! 
        {
            int k = 31; 
            if (argc > 2)
                k = stoi(argv[2]);
            cout << "building graph with k=" + to_string(k) + " for " + string(argv[1]) << endl;
            string args = "-in " + string(argv[1]) + " -kmer-size " + to_string(k) + " -abundance-min 1  -verbose 0  -max-memory 500 -mphf emphf";
            Graph graph = Graph::create (args.c_str());
            cout << "graph built, benchmarking.." << endl;
            debruijn_minim_bench (graph, k);
        }


    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }

    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }



}
