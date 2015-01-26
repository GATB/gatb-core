//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#undef NDEBUG
#include <cassert>

/********************************************************************************/
/*                          Fake nodes management                               */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We create an empty graph with a given kmer size
        Graph graph = Graph::create (7);

        // We create a sequence
        const char* seq = "ATGCATCGTA";

        // We ask for a fake node, starting at position 0 in the data
        Node n0 = graph.buildNode (seq+0);
        assert (graph.toString(n0) == "ATGCATC");

        // We ask for another fake node, starting at position 1 in the data
        Node n1 = graph.buildNode (seq+1);
        assert (graph.toString(n1) == "TGCATCG");

        std::cout << "Test OK" << std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
