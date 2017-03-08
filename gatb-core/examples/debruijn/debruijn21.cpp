//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*                   Looking for simple cycles in the graph.                    */
/*                                                                              */
/* Cmd-line: debruijn21 -graph <h5 file>                                        */
/*                                                                              */
/* Sample: debruijn21 -graph gatb-core/gatb-core/test/db/celegans_reads.h5      */
/*                                                                              */
/* Note:                                                                        */
/*     - '.h5' file contains the HDF5 formatted representation of a de bruijn   */
/*     graph created from a set of reads.                                       */
/*     - a '.h5' file is created using dbgh5 program provided with GATB-Core.   */
/*     Basic use is as follows:                                                 */
/*        dbgh5 -in <fasta/q file> -out <h5 file>                               */
/*     You can also control kmer-size and kmer abundance, see dbgh5 help.       */
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
