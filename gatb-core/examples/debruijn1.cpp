//! [snippet1]
// We include what we need for the test

#include <gatb/debruijn/impl/GraphFactory.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/math/NativeInt64.hpp>

// We use the required packages

using namespace std;

using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/

template<typename T>  void test (Graph<T>& graph)
{
    /** We get an iterator over all the nodes of the graph. */
    NodeIterator<T> itNodes = graph.nodes();

    /** We retrieve the first node. */
    itNodes.first();
    Node<T> node =  itNodes.item();

    /** We memorize the initial strand. */
    Strand strandDisplay = node.strand;

    /** We need an edges set. */
    EdgeSet<T> edges (graph);

    /** We loop until we reach a dead end along the kmers chain. */
    for (size_t nbEdges = ~0 ; nbEdges > 0; node = edges[0].to)
    {
        /** We retrieve the outcoming edges from the current node. */
        nbEdges = graph.getOutEdges (node, edges);

        /** If the current edge has some successors, we dump information about the first one.
         *  Note that we force the display in order to have the same strand for the two nodes.
         */
        if (nbEdges > 0)  {  cout << "-------------- " << edges[0].toString (strandDisplay) << " --------------" << endl; }

        /** The 'for' loop is going to update the current node with the first neighbor (if any). */
    }
}

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "you must provide:" << endl;
        cerr << "   1) kmer size"           << endl;
        cerr << "   2) solid kmers file"    << endl;
        cerr << "   3) debloom file"        << endl;
        return EXIT_FAILURE;
    }

    /** We create the graph. */
    Graph<NativeInt64> graph = GraphFactory::createGraph <NativeInt64> (
        new Property (STR_KMER_SIZE,   argv[1]),
        new Property (STR_KMER_SOLID,  argv[2]),
        new Property (STR_KMER_CFP,    argv[3]),
        PROP_END
    );

    /** We launch the test. */
    test<NativeInt64> (graph);
}
//! [snippet1]
