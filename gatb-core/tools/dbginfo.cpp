/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.push_back (new OptionOneParam (STR_URI_INPUT,  "graph file", true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We load the graph from the provided uri. */
        Graph graph = Graph::load (options->getStr(STR_URI_INPUT));

        /** We get information about the graph. */
        cout << graph.getInfo();
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
