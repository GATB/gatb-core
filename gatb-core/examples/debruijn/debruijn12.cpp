//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#undef NDEBUG
#include <cassert>

/********************************************************************************/
/*                          Node strands management                             */
/*                                                                              */
/* Cmd-line: debruijn12 (takes no argument)                                     */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We create the graph with a bank holding one sequence, and use a specific kmer size and solid kmer abundance to 1
        Graph graph = Graph::create (new BankStrings ("AATGC", NULL), "-kmer-size 5  -abundance-min 1  -verbose 0");

        // We get an iterator for all nodes of the graph.
        GraphIterator<Node> it = graph.iterator ();

        // We check that we have only one possible node
        assert (it.size() == 1);

        // We get the first (and only) node.
        it.first();  Node& current = it.item();

        // We get the ascii representation of the current iterated node
        std::string s = graph.toString (current);

        // We check the string value
        assert (s == "AATGC");

        // We reverse the node, which will be change its strand
        Node other = graph.reverse (current);

        // We check the string value of the reverse of the current node
        assert (graph.toString(other) == "GCATT");

        // We also check that the two nodes share the same kmer value
        assert (current.kmer == other.kmer);

        std::cout << "Test OK" << std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
