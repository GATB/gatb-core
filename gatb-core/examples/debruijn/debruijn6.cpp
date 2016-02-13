//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                      Graph nodes iteration                                   */
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

        // We get an iterator for all nodes of the graph.
        Graph::Iterator<Node> it = graph.iterator ();

        // We loop each node. Note the structure of the for loop.
        for (it.first(); !it.isDone(); it.next())
        {
            // The currently iterated node is available with it.item()

            // We dump an ascii representation of the current node.
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
