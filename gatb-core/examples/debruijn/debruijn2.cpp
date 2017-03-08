//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*  Graph creation from command line options                                    */
/*                                                                              */
/*  This snippet DO NOT use the OptionsParser facility of GATB-Core library     */
/*  but a simple string (expected to be a Fasta/q file) passed in on the        */
/*  command-line.                                                               */
/*                                                                              */
/* Cmd-line: debruijn1 <fasta/q file>                                           */
/*                                                                              */
/* Sample: debruijn1 gatb-core/gatb-core/test/db/reads1.fa                      */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be a bank URI).
    if (argc < 2)
    {
        std::cerr << "You must provide a bank uri." << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // We create the graph with the provided options.
        Graph graph = Graph::create ("-in %s", argv[1]);

        // We dump some information about the graph.
        std::cout << graph.getInfo() << std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
