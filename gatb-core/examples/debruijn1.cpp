//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/

void test (const Graph& graph)
{
    /** We get an iterator over all the nodes of the graph. */
    Graph::Iterator<Node> itNodes = graph.iterator<Node>();

    /** We retrieve the first node. */
    itNodes.first();
    Node node =  itNodes.item();

    /** We memorize the initial strand. */
    Strand strandDisplay = node.strand;

    /** We loop until we reach a dead end along the kmers chain. */
    for (Graph::Vector<Edge> edges; (edges = graph.successors<Edge>(node)).size() > 0; node = edges[0].to)
    {
        cout << "node=" << graph.toString(node, strandDisplay) << endl;
    }
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.add (new OptionOneParam (STR_URI_INPUT,  "graph file",               true));
    parser.add (new OptionOneParam (STR_KMER_SIZE,  "kmer size",                false, "27"));
    parser.add (new OptionOneParam (STR_NKS,        "kmer abundance threshold", false, "3" ));
    parser.add (new OptionNoParam  (STR_VERBOSE,    "verbosity",                false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

         /** We create the graph with the provided options. */
         Graph graph = Graph::create (options);

         /** We launch the test. */
         test (graph);

         /** We dump some information about the graph. */
         if (options->get(STR_VERBOSE) != 0)  {  std::cout << graph.getInfo() << endl;  }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
