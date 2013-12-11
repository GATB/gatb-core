//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#undef NDEBUG
#include <cassert>

/********************************************************************************/
/*                          Simple path iteration (with Edge)                   */
/********************************************************************************/
int main (int argc, char* argv[])
{
    size_t kmerSize = 11;
    char* seq = (char*) "AGGCGCTAGGGTAGAGGATGATGA";

    std::cout << "The initial sequence is '" << seq << "', kmer size is " << kmerSize << std::endl;

    // We create the graph from a given sequence, and for a given kmer size
    Graph graph = Graph::create (new BankStrings (seq, NULL),  "-kmer-size %d -nks 1", kmerSize);

    // We get the first node of the sequence.
    Node node = graph.buildNode (seq);

    // We create a Edge iterator that iterates all the simple node from the first node
    // Recall that a simple node has indegree==1 and outdegree==1
    Graph::Iterator<Edge> path = graph.simplePath<Edge> (node, DIR_OUTCOMING);

    // We iterate the simple path.
    for (path.first(); !path.isDone(); path.next())
    {
        std::cout << "   [" << path.rank() << "]  current item is " << graph.toString (path.item()) << std::endl;

        // We check that the transition nucleotide matches the correct character in the sequence
        assert (ascii(path.item().nt) == seq [graph.getKmerSize() + path.rank()]);
    }

    std::cout << "The simple path was " << path.rank() << " long" << std::endl;

    // We check that we found the correct number of nodes during the simple path iteration
    assert (path.rank() == strlen(seq) - graph.getKmerSize());

    return EXIT_SUCCESS;
}
//! [snippet1]
