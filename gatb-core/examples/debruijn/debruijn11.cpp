//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#undef NDEBUG
#include <cassert>

/********************************************************************************/
/*                          Node neighbors management                           */
/*                                                                              */
/* Cmd-line: debruijn11 (takes no argument)                                     */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We create the graph with a bank holding one sequence, and use a specific kmer size and solid kmer abundance to 1
        Graph graph = Graph::create (new BankStrings ("AATGC", NULL), "-kmer-size 4  -abundance-min 1  -verbose 0");

        // We get an iterator for all nodes of the graph.
        GraphIterator<Node> it = graph.iterator ();

        // We check that we have only two possible nodes
        assert (it.size() == 2);

        // We loop each node.
        for (it.first(); !it.isDone(); it.next())
        {
            // A shortcut.
            Node& current = it.item();

            // We get the ascii representation of the current iterated node
            std::string s = graph.toString (current);

            if (s=="AATG")
            {
                // We suppose that we know the only possible successor transition from this neighbor (nucl C)
                // It could have been retrieved by an previous call to Graph::neighbors<Edge> method

                // Now, we want to get the successor of the current node, only by giving the transition nucleotide.
                Node neighbor = graph.successor (current, NUCL_C);

                // WARNING ! This Graph::successor method doesn't check whether the neighbor actually belongs to the graph.
                // It is supposed here that the client knows perfectly that its transition nucleotide is valid.
                // This possibility is available for performance matters since checking graph membership may take some time.

                // We check the result
                assert (graph.toString (neighbor) == "ATGC");

                // There is an overloaded form for the Graph::successor method, where the user can provide an extra boolean.
                // With this implementation, a check is done to be sure whether the neighbor surely belongs to the graph.
                // If not so, the extra boolean is set to false, true otherwise.
                bool exists;
                Node potentialNeighbor;

                potentialNeighbor = graph.successor (current, NUCL_A, exists);
                assert (exists == false);

                potentialNeighbor = graph.successor (current, NUCL_C, exists);
                assert (exists == true);
                assert (graph.toString (potentialNeighbor) == "ATGC");

                potentialNeighbor = graph.successor (current, NUCL_G, exists);
                assert (exists == false);

                potentialNeighbor = graph.successor (current, NUCL_T, exists);
                assert (exists == false);
            }
        }

        std::cout << "Test OK" << std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
