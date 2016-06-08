//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                   Getting abundances for nodes graph.                        */
/*                                                                              */
/* This snippet shows how to retrieve the abundance of a node in the graph.     */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{

    cout << "This snippet shows how to retrieve the abundance of a node in the graph." << endl;

    if (argc < 3)
    {
        cout << "Arguments: [graphFile] [sequence]" << endl;
        exit(1);
    }

    char* graphFile = argv[1];
    char* seq       = argv[2];

    // We load the graph. IMPORTANT : must be have created with mphf (default parameter of dbgh5 nowadays)
    Graph graph = Graph::load (graphFile);

    // We build a fake node (we must be sure that it is in the graph).
    Node node = graph.buildNode (seq);

    // We query its abundance.
    cout << "abundance=" << graph.queryAbundance(node) << endl;
}
//! [snippet1]
