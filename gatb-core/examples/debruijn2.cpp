//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/

template<typename T>  T computeChecksum (Graph<T>& graph)
{
    /** We want to compute some checkum. */
    T checksum = 0;

    /** We get an iterator over all the nodes of the graph. */
    NodeIterator<T> itNodes = graph.nodes();

    NodeSet<T> nodes (graph);

    size_t nbNodesTotal      = 0;
    size_t nbSuccessorsTotal = 0;

    /** We iterate all the nodes of the graph. */
    for (itNodes.first(); !itNodes.isDone(); itNodes.next(), nbNodesTotal++)
    {
        /** We retrieve the successors. */
        size_t nbSuccessors = graph.getSuccessors (*itNodes, nodes);

        /** We update the number of found successors. */
        nbSuccessorsTotal += nbSuccessors;

        /** We iterate all the successors. */
        for (size_t i=0; i<nbSuccessors; i++)   {  checksum += nodes[i].kmer;  }
    }

    cout << "nbNodes=" << nbNodesTotal << "  nbSuccessors=" << nbSuccessorsTotal << "  checksum=" << checksum <<  endl;
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide:" << endl;
        cerr << "   1) reads file"  << endl;
        cerr << "   2) kmer size"   << endl;
        cerr << "optional:" << endl;
        cerr << "   3) nks"  << endl;
        return EXIT_FAILURE;
    }

    char*  bankUri  = argv[1];
    size_t kmerSize = atoi (argv[2]);
    size_t nks      = argc >=4 ? atoi (argv[3]) : 1;

    /** We create the graph with  1) a FASTA bank   2) a kmer size */
    Graph<NativeInt64> graph = GraphFactory::createGraph <NativeInt64> (new Bank (bankUri), kmerSize, nks);

    /** We want some time statistics. */
    TimeInfo ti;
    {
        TIME_INFO (ti, "checksum");

        /** We launch the test. */
        computeChecksum <NativeInt64> (graph);
    }

	/** We dump some information about the graph. */
    cout << graph.getInfo() << endl; 
    
    cout << "time=" << ti.getEntryByKey("checksum") << endl;
}
//! [snippet1]
