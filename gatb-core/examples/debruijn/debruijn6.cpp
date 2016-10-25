//! [snippet1]

// We include GATB-Core
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*            De Bruijn Graph nodes iteration                                   */
/********************************************************************************/
int main (int argc, char* argv[])
{
  // We check that the user provides at least one option: a Fasta/FastQ file.
  // Online GATB-Tutorial: this argument is automatically filled in with an 
  // appropriate file.
  if (argc < 2)
  {
    std::cerr << "Please, provide a sequence file." << std::endl;
    return EXIT_FAILURE;
  }

  try
  {
    // We load the graph from the given sequence file. In this example, we
    // create the graph only with k-mers observed at least 5 times in the data
    // set. We dot that using parameter '-abundance-min'.
    Graph graph = Graph::create (Bank::open(argv[1]), "-abundance-min %d", 5);

    // We get an iterator for all nodes of the graph.
    Graph::Iterator<Node> it = graph.iterator ();

    // We loop over each node. 
    for (it.first(); !it.isDone(); it.next())
    {
      // The currently iterated node is available with it.item().
      // Here, we use it just to dump an ascii representation of each node.
      std::cout << graph.toString (it.item()) << std::endl;
    }
  }
  catch (Exception& e)
  {
    std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
  }

  return EXIT_SUCCESS;
}
//! [snippet1]
