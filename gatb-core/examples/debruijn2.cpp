//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/

void computeChecksum (const Graph& graph, size_t nbCores)
{
    /** We want to compute some checkum. */
    struct Stats {
        Stats () : nbNodesTotal(0), nbSuccessorsTotal(0), nbBranching(0) {}
        size_t  nbNodesTotal;
        size_t  nbSuccessorsTotal;
        Integer checksum;
        size_t  nbBranching;
        void operator+= (const Stats& st)  {  checksum += st.checksum;   nbNodesTotal += st.nbNodesTotal;   nbSuccessorsTotal += st.nbSuccessorsTotal; }
    };

    ThreadObject<Stats> stats;

    /** We get an iterator over all the nodes of the graph. */
    Graph::Iterator<Node> itNodes = graph.iterator<Node>();

    /** We iterate all the nodes of the graph. */
    IDispatcher::Status status = Dispatcher(nbCores).iterate (itNodes, [&graph, &stats] (const Node& node)
    {
        Stats& s = stats();

        s.nbNodesTotal ++;

        /** We retrieve the successors. */
        Graph::Vector<Node> nodeset = graph.successors<Node> (node);

        /** We update the number of found successors. */
        s.nbSuccessorsTotal += nodeset.size();

        /** We iterate all the successors. */
        for (size_t i=0; i<nodeset.size(); i++)   {  s.checksum += nodeset[i].kmer;  }
    });

    stats.foreach ([&] (const Stats& st)  {  *stats += st;  });

    cout << "  nbNodes="      << stats->nbNodesTotal
         << "  nbSuccessors=" << stats->nbSuccessorsTotal
         << "  checksum="     << stats->checksum
         << "  time="         << status.time
         << "  nbCores="      << status.nbCores
         <<  endl;
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    const char* STR_NB_CORES_TEST = "-nb-cores-test";

    /** We create a command line parser. */
    OptionsParser parser;
    parser.add (new OptionOneParam (STR_URI_INPUT,     "graph file",               true));
    parser.add (new OptionOneParam (STR_KMER_SIZE,     "kmer size",                false, "27"));
    parser.add (new OptionOneParam (STR_NKS,           "kmer abundance threshold", false, "3" ));
    parser.add (new OptionOneParam (STR_NB_CORES_TEST, "nb cores (0 for all)",     false, "0"));
    parser.add (new OptionNoParam  (STR_VERBOSE,       "verbosity",                false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

         /** We create the graph with the provided options. */
         Graph graph = Graph::create (options);

         /** We launch the test. */
         for (size_t i=1; i<= options->getInt(STR_NB_CORES_TEST); i++)
         {
             computeChecksum (graph, i);
         }

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
