/* apart from te agreement part, it's becoming less and less useful given bench_mphf */

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


/* inspired by debruijn_test3 from unit tests*/

struct Parameter
{
    Parameter (size_t k, const Graph& graph) : graph(graph), k(k) {}
    const Graph& graph;
    size_t k;    
};


template<size_t span> struct debruijn_minim_bench {  void operator ()  (Parameter params)
{
    size_t kmerSize = params.k;
    const Graph& graph = params.graph;

    int miniSize = 8;
    int NB_REPETITIONS = 2000000;

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(3);

    GraphIterator<Node> nodes = graph.iterator ();

    cout << "graph has " << nodes.size() << " nodes" << endl;

    if (nodes.size() == 0)
        exit(1);

    nodes.first ();

    /** We get the first node. */
    Node node = nodes.item();

    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::ModelCanonical  ModelCanonical;
    typedef typename Kmer<span>::ModelDirect     ModelDirect;
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   ModelMini;
    typedef typename ModelMini::Kmer                        KmerType;

    ModelMini  modelMini (kmerSize, miniSize);
    
    /** We get the value of the current minimizer. */
    
    Type kmer = node.kmer.get<Type>();
    
    auto start_t=chrono::system_clock::now();
    for (unsigned int i = 0 ; i < NB_REPETITIONS ; i++)
        modelMini.getMinimizerValue(kmer, false);
    auto end_t=chrono::system_clock::now();
			
    cout << NB_REPETITIONS << " minimizers of length " << miniSize << " on a " << kmerSize << "-mer : " << diff_wtime(start_t, end_t) / unit << " seconds" << endl;
	
    cout << "---- now on all nodes of the graph -----\n";
    int times = max((int)(NB_REPETITIONS / nodes.size()), 1);


    /* compute a baseline */

    start_t=chrono::system_clock::now();
     for (int i=0; i < times; i++)
        for (nodes.first(); !nodes.isDone(); nodes.next())
            modelMini.getMinimizerValueDummy(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();
    auto baseline_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead (" << nodes.size() << " nodes, " << times << " times) : " << baseline_time << " seconds" << endl;

    /* existing code */

    start_t=chrono::system_clock::now();
    for (int i=0; i < times; i++)
        for (nodes.first(); !nodes.isDone(); nodes.next())
            modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), false);
    end_t=chrono::system_clock::now();
    cout << nodes.size() << " minimizers of length " << miniSize << " on all nodes (" << kmerSize << "-mers), " << to_string(times) << " times, with existing code : " << (diff_wtime(start_t, end_t) / unit) - baseline_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (int i=0; i < times; i++)
        for (nodes.first(); !nodes.isDone(); nodes.next())
            modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), true);
    end_t=chrono::system_clock::now();
    cout << nodes.size() << " minimizers of length " << miniSize << " on all nodes (" << kmerSize << "-mers), " << to_string(times) << " times, with new method    : " << (diff_wtime(start_t, end_t) / unit) - baseline_time << " seconds" << endl;
    //cout << modelMini._invalidMinimizersCounter << "/" << modelMini._minimizersCounter << " normal/fast minimizer computations" << endl,

   /* checking agreement between old and new method*/
    
    cout << "checking agreement... ";
    for (nodes.first(); !nodes.isDone(); nodes.next())
    {
        if (modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), false) != modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), true))
        {
            cout << "FAIL! problem with kmer " << graph.toString(nodes.item()) << " : (old) " << modelMini.getMinimizerString(nodes.item().kmer.get<Type>(), false) << " vs (new) " << modelMini.getMinimizerString(nodes.item().kmer.get<Type>(), true) << endl;
            cout << "debug: integer representation of new minimizer: " << modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), true) << endl;

            exit(1);
        }

        if (modelMini.getMinimizerPosition(nodes.item().kmer.get<Type>(), false) != modelMini.getMinimizerPosition(nodes.item().kmer.get<Type>(), true))
        {
            cout << "FAIL! problem with minimizer positions of kmer " << graph.toString(nodes.item()) << " : (old position) " << modelMini.getMinimizerPosition(nodes.item().kmer.get<Type>(), false) << " vs (new position) " << modelMini.getMinimizerPosition(nodes.item().kmer.get<Type>(), true) << endl;
            cout << " minimizer: " << modelMini.getMinimizerString(nodes.item().kmer.get<Type>()) << endl;

            exit(1);
        }

    }
    cout << "all good." << endl;
}
}; // end functor debruijn_minim_bench

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

            Integer::apply<debruijn_minim_bench, Parameter> (kmerSizes[j], Parameter( kmerSizes[j], graph) );

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
            Integer::apply<debruijn_minim_bench,Parameter> (k, Parameter(k, graph));
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
