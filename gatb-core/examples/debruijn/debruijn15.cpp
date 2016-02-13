//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                          Simple path iteration (with Edge)                   */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        size_t kmerSize = 11;
        char* seq = (char*) "AGGCGCTAGGGTAGAGGATGATGA";

        std::cout << "The initial sequence is '" << seq << "', kmer size is " << kmerSize << std::endl;

        // We create the graph from a given sequence, and for a given kmer size
        Graph graph = Graph::create (new BankStrings (seq, NULL),  "-kmer-size %d  -abundance-min 1  -verbose 0", kmerSize);

        // We get the first node of the sequence.
        Node node = graph.buildNode (seq);

        // We create a Edge iterator that iterates all the simple edges from the first node
        // Recall that a simple node has indegree==1 and outdegree==1
        Graph::Iterator<Edge> path = graph.simplePathEdge (node, DIR_OUTCOMING);

        // We iterate the simple path.
        for (path.first(); !path.isDone(); path.next())
        {
            std::cout << "   [" << path.rank() << "]  current item is " << graph.toString (path.item()) << std::endl;
        }

        std::cout << "The simple path was " << path.rank() << " long" << std::endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }


    return EXIT_SUCCESS;
}
//! [snippet1]
