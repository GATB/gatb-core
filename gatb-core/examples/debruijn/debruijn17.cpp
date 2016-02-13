//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                          Branching nodes neighbors                           */
/*                                                                              */
/* This snippet shows how to detect a specific pattern in the graph.            */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    size_t kmerSize = 7;

    // We define some sequences used for building our test graph.
    // Note that the sequences have a difference at index==kmerSize,
    // so the initial node AGGCGCT has 3 outcoming neighbors (and
    // therefore is a branching node)
    const char* sequences[] =
    {
        //      x <- difference here
        "AGGCGCTAGGGAGAGGATGATGAAA",
        "AGGCGCTCGGGAGAGGATGATGAAA",
        "AGGCGCTTGGGAGAGGATGATGAAA"
    };

    try
    {
        // We create the graph.
        Graph graph = Graph::create (new BankStrings (sequences, ARRAY_SIZE(sequences)),  "-kmer-size %d  -abundance-min 1  -verbose 0", kmerSize);

        // We get the first node (should be AGGCGCT); this is a branching node.
        Node node = graph.buildNode (sequences[0]);

        // We retrieve the branching neighbors for the node.
        Graph::Vector<BranchingEdge> branchingNeighbors = graph.successorsBranchingEdge (node);

        std::cout << "We found " << branchingNeighbors.size() << " branching neighbors from node " << graph.toString(node) << std::endl;

        // We loop over the branching neighbors. Here, we should have 3 branching neighbors, being the same GGGAGAG
        for (size_t i=0; i<branchingNeighbors.size(); i++)
        {
            // Note: we don't display all the transition nucleotides, only the first transition nucleotide.
            // We also display the number of transitions needed to link the two branching nodes.
            std::cout << graph.toString (branchingNeighbors[i])  << std::endl;
        }
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
