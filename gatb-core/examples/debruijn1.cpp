//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*             Graph creation from command line options                         */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We get a command line parser for graphs available options.
    OptionsParser parser = Graph::getOptionsParser();

    // We use a try/catch block in case we have some command line parsing issue.
    try  {
        // We parse the user options.
        parser.parse (argc, argv);

        // We create the graph with the provided options.
        Graph graph = Graph::create (parser.getProperties());

        // We dump some information about the graph.
        std::cout << graph.getInfo() << std::endl;
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
