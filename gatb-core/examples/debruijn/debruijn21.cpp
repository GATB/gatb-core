//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*                   Looking for simple cycles in the graph.                    */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be in HDF5 format).
    if (argc < 2)
    {
        cerr << "You must provide a HDF5 file." << endl;
        return EXIT_FAILURE;
    }

    // We create the graph with the bank and other options
    Graph graph = Graph::load (argv[1]);

    // We iterate the branching nodes
    Dispatcher().iterate (graph.iterator<BranchingNode> (), [&] (const BranchingNode& node)
    {
        // We iterate the successors of the current node
        graph.successors<BranchingEdge>(node).iterate ([&] (const BranchingEdge& edge)
        {
            if (edge.from == edge.to)  {  cout << "CYCLE: " << graph.toString (edge) << endl;  }
        });
    });

    return EXIT_SUCCESS;
}
//! [snippet1]
