//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*                  Path traversal in the graph                                 */
/*                                                                              */
/*  This snippet shows how to retrieve path in the graph with the Traversal     */
/*  class.                                                                      */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a command line parser.
    OptionsParser parser ("Traversal");
    parser.push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",  true));

    // We use a try/catch block in case we have some command line parsing issue.
    try
    {
        // We parse the user options.
        IProperties* props = parser.parse (argc, argv);

        // We load the graph
        Graph graph = Graph::load (props->getStr(STR_URI_GRAPH));

        // We create a Terminator object
        BranchingTerminator terminator (graph);

        // We create a Traversal instance for unitigs.
        Traversal* traversal = Traversal::create (TRAVERSAL_UNITIG, graph, terminator);
        LOCAL (traversal);

        Path pathRight;
        Path pathLeft;

        // We get one branching node in the graph
        Graph::Iterator<BranchingNode> itBranching = graph.iteratorBranching ();
        for (itBranching.first(); !itBranching.isDone(); itBranching.next())
        {
            BranchingNode current = itBranching.item();
            bool isMarked = terminator.is_marked_branching (current);

            cout << "Starting node : " << graph.toString (current) << "  marked=" << isMarked << endl;

            //if (terminator.is_marked_branching (itBranching.item()))  { continue; }

            int lenRight = traversal->traverse (current,                DIR_OUTCOMING, pathRight);
            Node rev = graph.reverse(current);
            int lenLeft  = traversal->traverse (rev, DIR_OUTCOMING, pathLeft);

            cout << "lenLeft=" << lenLeft << "  lenRight=" << lenRight << endl;
        }

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
