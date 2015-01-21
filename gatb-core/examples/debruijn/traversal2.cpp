//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

// We define a helper function to get the traversal kind from the command line
TraversalKind getTraversalKind (int argc, char* argv[]);

/********************************************************************************/
/*                  Path traversal in the graph                                 */
/*                                                                              */
/*  This snippet shows how to retrieve path in the graph with the Traversal     */
/*  class.                                                                      */
/*                                                                              */
/*  If one uses the 'unitig' mode, one should get only the first part before    */
/*  the SNP.                                                                    */
/*                                                                              */
/*  If one uses the 'contig' mode, one should get the whole sequence (which     */
/*  means a consensus sequence)                                                 */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We retrieve the traversal kind from the command line
    TraversalKind traversalKind = getTraversalKind (argc, argv);

//! [snippet1_traversal]
    const char* seqs[] =
    {
        "CGCTACAGCAGCTAGTTCATCATTGTTTATCAATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA",
        "CGCTACAGCAGCTAGTTCATCATTGTTTATCGATGATAAAATATAATAAGCTAAAAGGAAACTATAAATA"
        //      SNP HERE at pos 31      x
    };

    // We create a fake bank with a SNP
    IBank* bank = new BankStrings (seqs, ARRAY_SIZE(seqs));

    // We load the graph
    Graph graph = Graph::create (bank, "-abundance-min 1  -kmer-size 15  -verbose 0");

    // We create a Terminator object
    BranchingTerminator terminator (graph);

    // We create a Traversal instance according to the chosen traversal kind
    Traversal* traversal = Traversal::create (traversalKind, graph, terminator);
    LOCAL (traversal);

    // We create a node from the start of the first sequence
    Node node = graph.buildNode (seqs[0]);

    Path path;
    int len = traversal->traverse (node, DIR_OUTCOMING, path);

    // We dump the length, the starting node and the path
    cout << "length=" << len << " " << graph.toString (node) << path << endl;
//! [snippet1_traversal]

    return EXIT_SUCCESS;
}

/********************************************************************************/
TraversalKind getTraversalKind (int argc, char* argv[])
{
    const char* STR_TRAVERSAL_MODE = "-traversal";

    TraversalKind result;

    // We create a command line parser.
    OptionsParser parser ("Traversal");
    parser.push_back (new OptionOneParam (STR_TRAVERSAL_MODE, "traversal mode ('unitig' or 'contig'",  true));

    // We retrieve the traversal kind.
    try
    {
        IProperties* props = parser.parse (argc, argv);

        parse (props->getStr(STR_TRAVERSAL_MODE), result);
    }
    catch (OptionFailure& e)
    {
        e.displayErrors (std::cout);
        exit (EXIT_FAILURE);
    }
    catch (Exception& e)
    {
        cout << e.getMessage() << endl;
        exit (EXIT_FAILURE);
    }

    return result;
}

//! [snippet1]
