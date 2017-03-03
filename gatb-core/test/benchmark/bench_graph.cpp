/* this is the most advanced benchmark i've coded (more recent than bench_mphf, and some benchmarks overlap) 
 * */

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
#include <gatb/bank/api/IBank.hpp>


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

// for some reason it didn't recall those definition (yet in Minia.cpp, it did)
template <size_t span>
using NodeFast = Node_t<typename gatb::core::kmer::impl::Kmer<span>::Type >;
template <size_t span>
using EdgeFast = Edge_t<NodeFast<span> >;
template <size_t span>
using GraphDataVariantFast = boost::variant<GraphData<span>>;

struct Parameter
{
    Parameter (size_t k, string args, string seq="") : k(k), args(args), seq(seq){}
    size_t k; 
    string args;
    string seq;
};


template<size_t span> struct debruijn_mphf_bench {  void operator ()  (Parameter params)
{
    typedef NodeFast<span> NodeFastT;
    typedef GraphTemplate<NodeFastT,EdgeFast<span>,GraphDataVariantFast<span>> GraphFast;

    size_t kmerSize = params.k;
    
    Graph graph; 
    GraphFast graphFast;
  
    if (params.seq == "") 
    {
        graph = Graph::create (params.args.c_str());
        graphFast = GraphFast::create (params.args.c_str());
    }
    else
    {
        graph = Graph::create (new BankStrings (params.seq.c_str(), 0), params.args.c_str());
        graphFast = GraphFast::create (new BankStrings (params.seq.c_str(), 0), params.args.c_str());

    }

    cout << "graph built, benchmarking.." << endl;

    
    int miniSize = 8;
    int NB_REPETITIONS = 2000000;

    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(3);

    Graph::Iterator<Node> nodes = graph.iterator();
    typename GraphFast::template Iterator<NodeFastT> nodesFast = graphFast.iterator();
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

     // for some reason.. if *compiled*, this code confuses makes later MPHF queries 3x slower. really? yes. try to replace "if (confuse_mphf)" by "if (confuse_mphf && 0)" and re-run, you will see.
    {
        bool confuse_mphf = false;
        if (confuse_mphf)
        {
            //Type b; b.setVal(0); 
            //modelCanonical.emphf_hasher(modelCanonical.adaptor(b)); 
            //typedef std::pair<u_int8_t const*, u_int8_t const*> byte_range_t;

            //int c = 0; 
            //byte_range_t brange( reinterpret_cast <u_int8_t const*> (&c), reinterpret_cast <u_int8_t const*> (&c) + 2 );
            //byte_range_t brange( (u_int8_t const*) 1,(u_int8_t const*)33);
            //auto hashes = modelCanonical.empfh_hasher(brange);
        }


        for (int i = 0; i < 0; i++)
        {
            auto start_tt=chrono::system_clock::now();
            for (nodes.first(); !nodes.isDone(); nodes.next())
                modelCanonical.EMPHFhash(nodes.item().kmer.get<Type>());
            auto end_tt=chrono::system_clock::now();
            cout << "time to do " << nodes.size() << " computing EMPHFhash of kmers on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_tt, end_tt) / unit) << " seconds" << endl;
        }
        // it's slow. i don't understand why. see above for the "confuse mphf" part
        //return; //FIXME
    }


    /** We get the value of the first node (just an example, it's not used later). */
    Type kmer = node.kmer.get<Type>();
    
    auto start_t=chrono::system_clock::now();
    auto end_t=chrono::system_clock::now();
			
   cout << "----\non all nodes of the graph\n-----\n";

    /* disable node state (because we don't want to pay the price for overhea of checking whether a node is deleted or not in contain() */
   std::cout<< "PAY ATTENTION: this neighbor() benchmark, in the Bloom flavor, is without performing a MPHF query for each found node" << std::endl; 

   graph.disableNodeState(); 
   graphFast.disableNodeState(); 

   /* compute baseline times (= overheads we're not interested in) */

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
    {}
    end_t=chrono::system_clock::now();
    auto baseline_graph_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration (" << nodes.size() << " nodes) : " << baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
    {}
    end_t=chrono::system_clock::now();
    auto baseline_graphfast_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph NodeFast enumeration (" << nodes.size() << " nodes) : " << baseline_graphfast_time << " seconds" << endl;


    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelMini.getMinimizerValueDummy(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();
    auto baseline_minim_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration and minimizer computation setup (" << nodes.size() << " nodes) : " << baseline_minim_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        nodes.item().getKmer<Type>();
    end_t=chrono::system_clock::now();
    auto baseline_hash_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration and hash computation setup (" << nodes.size() << " nodes) : " << baseline_hash_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
       nodesFast.item().kmer;
    end_t=chrono::system_clock::now();
    auto baseline_hashfast_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph NodeFast enumeration and hash computation setup (" << nodes.size() << " nodes) : " << baseline_hashfast_time << " seconds" << endl;


    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.nodeMPHFIndexDummy(nodes.item());
    end_t=chrono::system_clock::now();
    auto baseline_mphf_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph nodes enumeration and mphf query setup (" << nodes.size() << " nodes) : " << baseline_mphf_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        graphFast.nodeMPHFIndexDummy(nodesFast.item());
    end_t=chrono::system_clock::now();
    auto baseline_mphffast_time = diff_wtime(start_t, end_t) / unit;
    cout << "baseline overhead for graph NodeFast enumeration and mphf query setup (" << nodes.size() << " nodes) : " << baseline_mphffast_time << " seconds" << endl;


    /* do actual benchmark */


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
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        graphFast.nodeMPHFIndex(nodesFast.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computations of MPHF index on all NodeFast (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_mphffast_time << " seconds" << endl;


    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelCanonical.getHash(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computing hash1 of kmers on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_hash_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelCanonical.getHash2(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computing hash2 of kmers on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_hash_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        modelCanonical.getHash2(nodesFast.item().kmer);
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computing hash2 of kmers on all NodeFast (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_hashfast_time << " seconds" << endl;


    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        modelCanonical.EMPHFhash(nodes.item().kmer.get<Type>());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " computing EMPHFhash of kmers on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_hash_time << " seconds" << endl;
    // it's slow. i don't understand why. see above for the "confuse mphf" part


    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.neighborsDummy(nodes.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " dummy neighbors() query on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_graph_time << " seconds" << endl;




    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.neighbors(nodes.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " neighbors() query on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        graphFast.neighbors(nodesFast.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " neighbors() query on all NodeFast (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_graphfast_time << " seconds" << endl;


/* isBranching */
    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.isBranching(nodes.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " isBranching() query on all nodes (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        graphFast.isBranching(nodesFast.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " isBranching() query on all NodeFast (" << kmerSize << "-mers) : " << (diff_wtime(start_t, end_t) / unit) - baseline_graphfast_time << " seconds" << endl;



    
    /* now, compute adjacency! */



    graph.precomputeAdjacency();
    graphFast.precomputeAdjacency();

    cout << "adjacency precomputed" << endl;

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.neighbors(nodes.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " neighbors() query on all nodes (" << kmerSize << "-mers) using adjacency : " << (diff_wtime(start_t, end_t) / unit) - baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        graphFast.neighbors(nodesFast.item());
    end_t=chrono::system_clock::now();
    cout << "time to do " << nodes.size() << " fast neighbors() query on all NodeFast (" << kmerSize << "-mers) using adjacency : " << (diff_wtime(start_t, end_t) / unit) - baseline_graphfast_time << " seconds" << endl;
    
    /* isBranching */

    start_t=chrono::system_clock::now();
    for (nodes.first(); !nodes.isDone(); nodes.next())
        graph.isBranching(nodes.item());
    end_t=chrono::system_clock::now();

    cout << "time to do " << nodes.size() << " isBranching() query on all nodes (" << kmerSize << "-mers) using adjacency : " << (diff_wtime(start_t, end_t) / unit) - baseline_graph_time << " seconds" << endl;

    start_t=chrono::system_clock::now();
    for (nodesFast.first(); !nodesFast.isDone(); nodesFast.next())
        graphFast.isBranching(nodesFast.item());
    end_t=chrono::system_clock::now();
    cout << "time to do " << nodes.size() << " fast isBranching() query on all NodeFast (" << kmerSize << "-mers) using adjacency : " << (diff_wtime(start_t, end_t) / unit) - baseline_graphfast_time << " seconds" << endl;
    


    /** We remove the graph. */
    //graph.remove ();
    //graphFast.remove (); // no actually, I want to keep the .h5 file


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
    size_t kmerSizes[] = {31};

    for (size_t i=0; i<ARRAY_SIZE(sequences); i++)
    {
        for (size_t j=0; j<ARRAY_SIZE(kmerSizes); j++)
        {
            string args = "-kmer-size " + std::to_string(kmerSizes[j]) + "  -abundance-min 1  -verbose 0  -max-memory 500";

            Integer::apply<debruijn_mphf_bench, Parameter> (kmerSizes[j], Parameter( kmerSizes[j], args, sequences[i]) );
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
            Integer::apply<debruijn_mphf_bench, Parameter> (k, Parameter( k, args) );
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
