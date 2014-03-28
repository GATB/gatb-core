//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*               */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be in HDF5 format).
    if (argc < 2)
    {
        cerr << "You must provide a HDF5 file." << endl;
        return EXIT_FAILURE;
    }

    // We create the graph with the bank and other options
    Graph graph = Graph::load (argv[1]);

    // We iterate the branching nodes
    graph.iterator<BranchingNode> ().iterate ([&] (BranchingNode& node)
    {
        // We get the branching neighbors of the current node
        Graph::Vector<BranchingEdge> neighbors = graph.successors<BranchingEdge> (node);

        // We look for
        for (size_t i=0; i<neighbors.size(); i++)
        {
            for (size_t j=i+1; j<neighbors.size(); j++)
            {
                if (neighbors[i].to == neighbors[j].to)
                {
                    cout << "BUBBLE" << endl;
                    {
                        cout << "   " << graph.toString (neighbors[i].from);
                        Nucleotide nt = neighbors[i].nt;
                        Node next = graph.successor<Node> (node, nt);

                        while (1)
                        {
                            cout << ascii (nt);
                            Graph::Vector<Edge> simpleNeighbors = graph.successors<Edge> (next);
                            if (simpleNeighbors.size() != 1)  { break; }
                            next = simpleNeighbors[0].to;
                            nt = simpleNeighbors[0].nt;
                        }
                        cout << endl;
                    }
                    {
                        cout << "   " << graph.toString (neighbors[j].from);
                        Nucleotide nt = neighbors[j].nt;
                        Node next = graph.successor<Node> (node, nt);

                        while (1)
                        {
                            cout << ascii (nt);
                            Graph::Vector<Edge> simpleNeighbors = graph.successors<Edge> (next);
                            if (simpleNeighbors.size() != 1)  { break; }
                            next = simpleNeighbors[0].to;
                            nt = simpleNeighbors[0].nt;
                        }
                        cout << endl;
                    }
                }
            }
        }
    });

    return EXIT_SUCCESS;
}
//! [snippet1]
