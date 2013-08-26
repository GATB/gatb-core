//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

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

    /** We take the first node of the graph. */
    Node<NativeInt64> node = graph.nodes().item();

    /** We get the branching interval for this node. */
    Node<NativeInt64> begin, end;   graph.getNearestBranchingRange (node, begin, end);

    /** We dump the interval nodes. */
    cout << begin.toString() << "  "  << end.toString()  << endl;
}
//! [snippet1]
