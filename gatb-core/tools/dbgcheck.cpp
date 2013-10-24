/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/

struct Stats
{
    Stats() :  nbSuccessors(0), abundance(0) {}
    Integer    checksumNodes;
    Integer    checksumSuccessors;
    u_int64_t  nbSuccessors;
    u_int64_t  abundance;

    Stats& operator+= (const Stats& st)
    {
        if (this != &st)
        {
            checksumNodes      += st.checksumNodes;
            checksumSuccessors += st.checksumSuccessors;
            nbSuccessors       += st.nbSuccessors;
            abundance          += st.abundance;
        }
        return *this;
    }
};

/********************************************************************************/

int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.add (new OptionOneParam (STR_URI_INPUT,  "graph file", true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We load the graph from the provided uri. */
        Graph graph  = GraphFactory::load (options->getStr(STR_URI_INPUT));

        /** We get an iterator over the nodes of the graph. */
        INodeIterator* nodes = graph.nodes ();

        /** We want to gather some statistics during the iteration. */
        ThreadContainer<Stats> stats;

        /** We launch the iteration.
         *  Note that we encapsulate the node iterator by an iterator that writes progression status on the console. */
        IDispatcher::Status status = ParallelDispatcher().iterate (ProgressIterator<ProgressTimer>(nodes), [&graph,&stats] (const Node& node)
        {
            /** We get the Stats object for the current thread. */
            Stats& s = stats.current();

            /** We update the statistics. */
            s.checksumNodes  += node.kmer.value;
            s.abundance      += node.kmer.abundance;

            /** We retrieve the successors. */
            NodeSet nodeset (graph);
            size_t nbSuccessors = graph.getSuccessors (node, nodeset);

            s.nbSuccessors += nbSuccessors;

            /** We iterate all the successors. */
            for (size_t i=0; i<nbSuccessors; i++)  {  s.checksumSuccessors += nodeset[i].kmer.value;  }
        });

        /** We aggregate the statistics gathered during iteration. */
        Stats s;   stats.foreach ([&s] (const Stats& st)  {  s += st;  });

        /** We dump the statistics. */
        std::cout << std::endl;
        std::cout << "nbSuccessors       = " << s.nbSuccessors           << "  "  << std::endl
                  << "checkumNodes       = " << s.checksumNodes          << "  "  << std::endl
                  << "checksumSuccessors = " << s.checksumSuccessors     << "  "  << std::endl
                  << "abundance          = " << s.abundance              << "  "  << std::endl
                  << "time               = " << status.time              << "  "  << std::endl
                  << "nbCores            = " << status.nbCores           << "  "  << std::endl
                  << std::endl;
        std::cout << std::endl;
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
