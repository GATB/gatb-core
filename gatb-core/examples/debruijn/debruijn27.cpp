//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                   Getting abundances for nodes graph.                        */
/*                                                                              */
/* This snippet shows how to retrieve the abundance of a node in the graph.     */
/*                                                                              */
/* Cmd-line: debruijn27 <h5 file> <sequence>                                    */
/*                                                                              */
/* Sample: debruijn27 gatb-core/gatb-core/test/db/celegans_reads.h5 ...         */
/*                    "TTTGCCCATTTCCTGCCATTTGTC"                                */
/*                     |-> 5' end of 1st sequence from                          */
/*                        https://github.com/GATB/gatb-core-tuto/blob/master/server/data/celegans_reads.fasta*/
/*                                                                              */
/* Note:                                                                        */
/*     - '.h5' file contains the HDF5 formatted representation of a de bruijn   */
/*     graph created from a set of reads.                                       */
/*     - a '.h5' file is created using dbgh5 program provided with GATB-Core.   */
/*     Basic use is as follows:                                                 */
/*        dbgh5 -in <fasta/q file> -out <h5 file>                               */
/*     You can also control kmer-size and kmer abundance, see dbgh5 help.       */
/********************************************************************************/
int main (int argc, char* argv[])
{

    cout << "This snippet shows how to retrieve the abundance of a node in the graph." << endl;

    if (argc < 3)
    {
        cout << "Arguments: [graphFile] [sequence]" << endl;
        exit(1);
    }

    char* graphFile = argv[1];
    char* seq       = argv[2];

    // We load the graph. IMPORTANT : must be have created with mphf (default parameter of dbgh5 nowadays)
    Graph graph = Graph::load (graphFile);

    // We build a fake node (we must be sure that it is in the graph).
    Node node = graph.buildNode (seq);

    // We query its abundance.
    cout << "abundance=" << graph.queryAbundance(node) << endl;
}
//! [snippet1]
