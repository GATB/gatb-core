//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*             Graph creation from a bank and command line options              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be a FASTA file URI).
    if (argc < 2)
    {
        std::cerr << "You must provide a FASTA file uri." << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // We create the graph with the bank and other options
        Graph graph = Graph::create (Bank::open(argv[1]), "-abundance-min %d", 5);

        // We dump some information about the graph.
        std::cout << graph.getInfo() << std::endl;

        // Note: Graph::create will take care about 'bank' object and will delete it if nobody else needs it.
        // In other words: there is no need here to call 'delete' on 'bank' here.
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
