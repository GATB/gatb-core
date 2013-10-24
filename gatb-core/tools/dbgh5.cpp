/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/gatb_core.hpp>

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.add (new OptionOneParam (STR_URI_INPUT,  "reads input",              true));
    parser.add (new OptionOneParam (STR_URI_OUTPUT, "graph output",             false));
    parser.add (new OptionOneParam (STR_KMER_SIZE,  "kmer size",                false, "27"));
    parser.add (new OptionOneParam (STR_NKS,        "kmer abundance threshold", false, "3" ));
    parser.add (new OptionOneParam (STR_MAX_MEMORY, "max memory",               false, "1000"));
    parser.add (new OptionOneParam (STR_MAX_DISK,   "max disk",                 false, "0"));
    parser.add (new OptionOneParam (STR_NB_CORES,   "nb cores (0 for all)",     false, "0"));
    parser.add (new OptionNoParam  (STR_VERBOSE,    "verbosity",                false));
    parser.add (new OptionNoParam  (STR_HELP,       "help",                     false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

         /** We create the graph with the provided options. */
         Graph graph = GraphFactory::create (options);

         /** We dump some information about the graph. */
         if (options->get(STR_VERBOSE) != 0)  {  std::cout << graph.getInfo() << std::endl;  }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
