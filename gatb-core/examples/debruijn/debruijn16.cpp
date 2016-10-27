//! [snippet1]
// GATB-Core Online Tutorial 

/********************************************************************************/
/*                          Branching nodes analysis                            */
/*                                                                              */
/* This snippet shows how to detect a specific pattern in the De Bruijn graph.  */
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>

int main (int argc, char* argv[])
{
  // We set the kmer size
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
    // We create the graph from our sample set of sequences.
    // Graph is created with the following parameters:
    //  kmer-size is 7 (see declaration, above)
    //  we retain kmers for which counting (i.e. abundance) is 1
    //  we request GATB-Core to be silent while creating the graph
    Graph graph = Graph::create (
      new BankStrings (sequences, ARRAY_SIZE(sequences)),  
      "-kmer-size %d  -abundance-min 1  -verbose 0", kmerSize);

    // We get the first node (should be AGGCGCT); this is a branching node.
    Node node = graph.buildNode (sequences[0]);

    // We retrieve the branching neighbours for that node.
    Graph::Vector<BranchingNode> branchingNeighbors = graph.successorsBranching (node);

    std::cout << "We found " << branchingNeighbors.size() << " branching neighbors from node " << graph.toString(node) << std::endl;

    // We loop over the branching neighbors. 
    // Here, we should have 3 branching neighbors, all being the same: GGGAGAG
    for (size_t i=0; i<branchingNeighbors.size(); i++)
    {
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
