//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*             Graph loading from a HDF5 file                                   */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be in HDF5 format).
    if (argc < 2)
    {
        std::cerr << "You must provide a HDF5 file." << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // We create the graph with the bank and other options
        Graph graph = Graph::load (argv[1]);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
