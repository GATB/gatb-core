//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/* Graph creation from a fake bank and no command line options                  */
/*                                                                              */
/* You can dump the result of this test with HDF5 tools (provided by GATB):     */
/*      h5dump mygraph.h5                                                       */
/*                                                                              */
/* Cmd-line: debruijn4 (takes no argument)                                      */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We get a handle on a fake bank made of 3 sequences.
    IBank* bank = new BankStrings (
        "ATCGTACGACGCTAGCTAGCA",
        "ACTACGTATCGGTATATATTTTCGATCGATCAG",
        "TGACGGTAGCATCGATCAGGATCGA",
        NULL
    );

    try
    {
        // We create the graph with the bank and other options
        Graph graph = Graph::create (bank, "-kmer-size 5  -abundance-min 1  -out mygraph");

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
