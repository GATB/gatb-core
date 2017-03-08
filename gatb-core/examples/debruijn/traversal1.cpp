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
/* Cmd-line: traversal1 -graph <h5 file>                                        */
/*                                                                              */
/* Sample: traversal1 -graph gatb-core/gatb-core/test/db/celegans_reads.h5      */
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
        GraphIterator<BranchingNode> itBranching = graph.iteratorBranching ();
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
