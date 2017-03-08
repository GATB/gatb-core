//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/* Graph loading from a HDF5 file, the iterate over the branching nodes.        */
/*                                                                              */
/* Cmd-line: debruijn8 <h5 file>                                                */
/*                                                                              */
/* Sample: debruijn8 gatb-core/gatb-core/test/db/celegans_reads.h5              */
/*                                                                              */
/* Note:                                                                        */
/*     - '.h5' file contains the HDF5 formatted representation of a de bruijn   */
/*     graph created from a set of reads.                                       */
/*     - a '.h5' file is created using dbgh5 program provided with GATB-Core.   */
/*     Basic use is as follows:                                                 */
/*        dbgh5 -in <fasta/q file> -out <h5 file>                               */
/*     You can also control kmer-size and kmer abundance, see dbgh5 help.       */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be in HDF5 format).
    if (argc < 2)
    {
        std::cerr << "You must provide a HDF5 file." << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // We load the graph from the given graph file
        Graph graph = Graph::load (argv[1]);

        // We get an iterator for branching nodes of the graph.
        GraphIterator<BranchingNode> it = graph.iteratorBranching ();

        // We loop each node. Note the structure of the for loop.
        for (it.first(); !it.isDone(); it.next())
        {
            // The currently iterated branching node is available with it.item()

            // We dump an ascii representation of the current node.
            std::cout << "[" << it.rank() << "] " << graph.toString (it.item()) << std::endl;
        }
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
