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
    /** We create a command line parser. */
    OptionsParser parser ("GraphStats");
    parser.push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",  true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We load the graph
        Graph graph = Graph::load (options->getStr(STR_URI_GRAPH));

        // We iterate the branching nodes
        Dispatcher().iterate (graph.iteratorBranching (), [&] (const BranchingNode& node)
        {
            // We iterate the successors of the current node
            graph.successorsBranchingEdge((Node&)node).iterate ([&] (const BranchingEdge& edge)
            {
                if (edge.from == edge.to)  {  cout << "CYCLE: " << graph.toString (edge) << endl;  }
            });
        });
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
