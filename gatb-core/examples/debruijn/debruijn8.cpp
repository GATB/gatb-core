//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                      Graph branching nodes iteration                         */
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
        Graph::Iterator<BranchingNode> it = graph.iteratorBranching ();

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
