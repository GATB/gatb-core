//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/

template<typename T>  void test (Graph<T>& graph)
{
    /** We get an iterator over all the nodes of the graph. */
    NodeIterator<T> itNodes = graph.nodes();

    /** We retrieve the first node. */
    itNodes.first();
    Node<T> node =  itNodes.item();

    /** We memorize the initial strand. */
    Strand strandDisplay = node.strand;

    /** We need an edges set. */
    EdgeSet<T> edges (graph);

    /** We loop until we reach a dead end along the kmers chain. */
    for (size_t nbEdges = ~0 ; nbEdges > 0; node = edges[0].to)
    {
        /** We retrieve the outcoming edges from the current node. */
        nbEdges = graph.getOutEdges (node, edges);

        /** If the current edge has some successors, we dump information about the first one.
         *  Note that we force the display in order to have the same strand for the two nodes.
         */
        if (nbEdges > 0)  {  cout << "-------------- " << edges[0].toString (strandDisplay) << " --------------" << endl; }

        /** The 'for' loop is going to update the current node with the first neighbor (if any). */
    }
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

    /** We launch the test. */
    test<NativeInt64> (graph);
}
//! [snippet1]
