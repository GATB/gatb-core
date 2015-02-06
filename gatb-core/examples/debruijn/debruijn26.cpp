//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                   Getting abundances for nodes graph.                        */
/*                                                                              */
/* This snippet shows how to retrieve the abundance of a node in the graph.     */
/* This feature is enabled only when the "-mphf emphf" is set during graph      */
/* creation.                                                                    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    const char* seq = "AAAAACTACATTACCCGTTTGCGAGACAGGTA";

    size_t kmerSize = strlen (seq) - 1;

    // We create a fake bank with 4 times the same sequence
    IBank* bank = new BankStrings (seq, seq, seq, seq, 0);

    // We create the graph. IMPORTANT : we use the -mphf option
    Graph graph = Graph::create (bank, "-kmer-size %d  -abundance-min 1  -mphf emphf  -verbose 0", kmerSize);

    // We build a fake node (we are sure that it will be in the graph).
    Node node = graph.buildNode (seq);

    // We query its abundance.
    cout << "abundance=" << graph.queryAbundance(node) << endl;
}
//! [snippet1]
