
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


template<size_t span> struct debruijn_mphf_bench {  void operator ()  (Parameter params)
{
    size_t kmerSize = params.k;
    const Graph& graph = params.graph;
    
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
    ModelCanonical  modelCanonical (kmerSize);

    /** We get the value of the current minimizer. */
    Type kmer = node.kmer.get<Type>();
    
    auto start_t=chrono::system_clock::now();
    for (unsigned int i = 0 ; i < NB_REPETITIONS ; i++)
        modelMini.getMinimizerValue(kmer);
    auto end_t=chrono::system_clock::now();
			
    cout << "time to do " << NB_REPETITIONS << " computations of minimizers of length " << miniSize << " on a " << kmerSize << "-mer : " << diff_wtime(start_t, end_t) / unit << " seconds" << endl;
	
    start_t=chrono::system_clock::now();
    for (unsigned int i = 0 ; i < NB_REPETITIONS ; i++)
         graph.nodeMPHFIndex(node);
    end_t=chrono::system_clock::now();
    //cout << "mphf index of node " << graph.toString(node) << " is " << graph.nodeMPHFIndex(node) << endl;

    cout << "time to do " << NB_REPETITIONS << " computations of MPHF index on a " << kmerSize << "-mer : " << diff_wtime(start_t, end_t) / unit << " seconds" << endl;
		
    cout << "----\nnow on all nodes of the graph\n-----\n";

    /* compute baseline times (= overheads we're not interested in) */

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
    {}
    end_t=chrono::system_clock::now();
    auto baseline_graph_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration (" << nodes.size() << " nodes) : " << baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelMini.getMinimizerValueDummy(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();
    auto baseline_minim_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration and minimizer computation setup (" << nodes.size() << " nodes) : " << baseline_minim_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.nodeMPHFIndexDummy(nodes.item());
    end_t=chrono::system_clock::now();
    auto baseline_mphf_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration and mphf query setup (" << nodes.size() << " nodes) : " << baseline_mphf_time << " seconds" << endl;

    /* do actual benchmark */


    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), false);
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computations of minimizers of length " << miniSize << " on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_minim_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelMini.getMinimizerValue(nodes.item().kmer.get<Type>(), true);
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computations of minimizers (fast method) of length " << miniSize << " on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_minim_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.nodeMPHFIndex(nodes.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computations of MPHF index on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_mphf_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelCanonical.reverse(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " revcomps of all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelMini.sweepForAA(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " just sweepings for AA's across kmers on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_minim_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelCanonical.getHash(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computing hash1 of kmers on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_minim_time << " seconds" << endl;



}
};

void debruijn_mphf ()
{
    const char* sequences [] =
    {
//        "ACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAAT",
//        "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA",
        "CGCTATTCATCATTGTTTATCAATGAGCTAAAAGGAAACTATAAATAACCATGTATAATTATAAGTAGGTACCTATTTTTTTATTTTAAACTGAAATTCAATATTATATAGGCAAAG"
    };

    //        size_t kmerSizes[] = { 13, 15, 17, 19, 21, 23, 25, 27, 29, 31};
    size_t kmerSizes[] = {81};

    for (size_t i=0; i<ARRAY_SIZE(sequences); i++)
    {
        for (size_t j=0; j<ARRAY_SIZE(kmerSizes); j++)
        {
            /** We create the graph. */
            Graph graph = Graph::create (
                    new BankStrings (sequences[i], 0),
                    "-kmer-size %d  -abundance-min 1  -verbose 0  -max-memory %d", kmerSizes[j], 500 
                    );

            Integer::apply<debruijn_mphf_bench, Parameter> (kmerSizes[j], Parameter( kmerSizes[j], graph) );

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
            debruijn_mphf();
        else
            // else, use a real file! 
        {
            cout << "building graph" << endl;
            int k = 31;
            if (argc > 2)
                k = stoi(argv[2]);

            string args = "-in " + string(argv[1]) + " -kmer-size " + std::to_string(k) + " -abundance-min 1  -verbose 0  -max-memory 500";
            Graph graph = Graph::create (args.c_str());
            cout << "graph built, benchmarking.." << endl;
            Integer::apply<debruijn_mphf_bench, Parameter> (k, Parameter( k, graph) );
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
