//! [snippet1]
// We include what we need for the test

#include <gatb/debruijn/impl/GraphFactory.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/math/NativeInt64.hpp>

// We use the required packages

using namespace std;

using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/

template<typename T>  T computeChecksum (Graph<T>& graph)
{
    /** We want to compute some checkum. */
    T checksum = 0;

    /** We get an iterator over all the nodes of the graph. */
    NodeIterator<T> itNodes = graph.nodes();

    NodeSet<T> nodes (graph);

    size_t nbNodesTotal      = 0;
    size_t nbSuccessorsTotal = 0;

    /** We iterate all the nodes of the graph. */
    for (itNodes.first(); !itNodes.isDone(); itNodes.next(), nbNodesTotal++)
    {
        /** We retrieve the successors. */
        size_t nbSuccessors = graph.getSuccessors (*itNodes, nodes);

        /** We update the number of found successors. */
        nbSuccessorsTotal += nbSuccessors;

        /** We iterate all the successors. */
        for (size_t i=0; i<nbSuccessors; i++)   {  checksum += nodes[i].kmer;  }
    }

    cout << "nbNodes=" << nbNodesTotal << "  nbSuccessors=" << nbSuccessorsTotal << "  checksum=" << checksum <<  endl;
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "you must provide:" << endl;
        cerr << "   1) kmer size"                    << endl;
        cerr << "   2) solid kmers file"             << endl;
        cerr << "   3) critical false positive file" << endl;
        return EXIT_FAILURE;
    }

    /** We create the graph with the user parameters. */
    Graph<NativeInt64> graph = GraphFactory::createGraph <NativeInt64> (
        new Property (STR_KMER_SIZE,   argv[1]),
        new Property (STR_KMER_SOLID,  argv[2]),
        new Property (STR_KMER_CFP,    argv[3]),
        PROP_END
    );

    /** We launch the test. */
    computeChecksum <NativeInt64> (graph);
}
//! [snippet1]
