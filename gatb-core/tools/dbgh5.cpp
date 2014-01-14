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
    OptionsParser parser = Graph::getOptionsParser();
    parser.setName ("dbgh5");

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

         /** We create the graph with the provided options. */
         Graph graph = Graph::create (options);

         /** We dump some information about the graph. */
         if (options->get(STR_VERBOSE) != 0)  {  std::cout << graph.getInfo() << std::endl;  }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        e.getParser().displayVersion(stdout);
        return EXIT_FAILURE;
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
    catch (char const* msg)
    {
        std::cout << "EXCEPTION: " << msg << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
