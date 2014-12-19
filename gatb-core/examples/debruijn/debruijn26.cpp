//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/********************************************************************************/
int main (int argc, char* argv[])
{
    const char* seq = "AAAAACTACATTACCCGTTTGCGAGACAGGTA";

    size_t kmerSize = strlen (seq) - 1;

    // We create a fake bank
    IBank* bank = new BankStrings (seq, seq, seq, seq, 0);

    // We create the graph. IMPORTANT : we use the -mphf option
    Graph graph = Graph::create (bank, "-kmer-size %d  -abundance-min 1  -mphf emphf", kmerSize);

    // We build a fake node (we are sure that it will be in the graph).
    Node node = graph.buildNode ((char*)seq, 0);

    // We query its abundance.
    cout << "abundance=" << graph.queryAbundance(node) << endl;
}
//! [snippet1]
