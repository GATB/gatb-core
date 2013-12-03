//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*             Graph creation from command line string                          */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be a FASTA file URI).
    if (argc < 2)
    {
        std::cerr << "You must provide a FASTA file uri." << std::endl;
        return EXIT_FAILURE;
    }

    // We create the graph with the provided options.
    Graph graph = Graph::create ("-in %s", argv[1]);

    // We dump some information about the graph.
    std::cout << graph.getInfo() << std::endl;

    return EXIT_SUCCESS;
}
//! [snippet1]
