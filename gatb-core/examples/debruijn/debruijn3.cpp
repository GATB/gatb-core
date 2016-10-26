//! [snippet1]
// GATB-Core Online Tutorial 
/********************************************************************************/
/*             Graph creation from a bank file                                  */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>

/********************************************************************************/
// START Application
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

    // We dump some information about the graph.
    std::cout << graph.getInfo() << std::endl;

    // Note: Graph::create will take care about 'bank' object and will delete it 
    // if nobody else needs it.
    // In other words: there is no need here to call 'delete' on 'bank' here.
    }
    catch (Exception& e)
    {
      std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }
    return EXIT_SUCCESS;
}
//! [snippet1]
