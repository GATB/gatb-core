//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <map>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <set>

#undef NDEBUG
#include <cassert>

using namespace std;

/********************************************************************************/

template<typename T>
class GraphMarker
{
public:

    GraphMarker (const Graph& graph) : graph(graph)
    {
        graph.iterator<T>().iterate ([&] (const T& item)  {  this->markMap[item] = false;  });
    }

    void mark (const T& item)  {  markMap [item] = true;  }

    bool isMarked (const T& item) const  {  return markMap.find(item)->second;  }

private:
    const Graph& graph;
    map<T,bool>  markMap;
};

/********************************************************************************/

template<typename T>
class BFS
{
public:

    BFS (const Graph& graph, GraphMarker<T>& marker, Direction dir) : graph(graph), marker(marker), direction(dir)  {}

    u_int64_t run (const T& node)
    {
        // We initialize the queue used for the recursion
        while (!frontline.empty())  { frontline.pop(); }
        frontline.push (node);

        // We reset the already seen nodes.
        alreadyFrontlined.clear();

        // We add the initial node to the set of seen nodes.
        alreadyFrontlined.insert(node.kmer);

        // We launch the BFS visit.
        while (nextDepth() == true)  {}

        // We return the number of nodes for this connected component
        return alreadyFrontlined.size();
    }

private:
    const Graph&     graph;
    GraphMarker<T>&  marker;
    Direction        direction;
    queue<T>         frontline;
    set<Node::Value> alreadyFrontlined;

    bool nextDepth ()
    {
        queue<T> newFrontline;

        while (!frontline.empty())
        {
            T current = frontline.front();
            frontline.pop();

            // We mark the node
            marker.mark (current);

            /** We loop the neighbors of the current node. */
            Graph::Vector<T> neighbors = graph.neighbors<T> (current, direction);

            for (size_t i=0; i<neighbors.size(); i++)
            {
                // if this bubble contains a marked (branching) kmer, stop everyone at once (to avoid redundancy)
                if (marker.isMarked (neighbors[i]))  {  continue;  }

                // test if that node hasn't already been explored
                if (alreadyFrontlined.find (neighbors[i].kmer) != alreadyFrontlined.end())  { continue; }

                /** We add the new node to the new front line. */
                newFrontline.push (neighbors[i]);

                /** We memorize the new node. */
                alreadyFrontlined.insert (neighbors[i].kmer);
            }
        }

        frontline = newFrontline;

        return frontline.empty() == false;
    }

};

/********************************************************************************/
/*        Computing connected components of the branching nodes subgraph.       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides at least one option (supposed to be in HDF5 format).
    // IMPORTANT: only the prefix of the file has to be given. For instance, if the graph
    // file name is "foo.h5", one should provide "foo"
    if (argc < 2)
    {
        cerr << "You must provide a HDF5 file." << endl;
        return EXIT_FAILURE;
    }

    // We create the graph with the bank and other options
    Graph graph = Graph::load (argv[1]);

    // We create a graph marker.
    GraphMarker<BranchingNode> marker (graph);

    // We create an object for Breadth First Search for the de Bruijn graph.
    BFS<BranchingNode> bfs (graph, marker, DIR_OUTCOMING);

    // We want to compute the distribution of connected components of the branching nodes.
    //    - key is a connected component class (for a given number of branching nodes for this component)
    //    - value is the number of times this component class occurs in the branching sub graph
    map<size_t,size_t> distrib;

    // We get an iterator for all nodes of the graph. We use a progress iterator to get some progress feedback
    ProgressIterator<BranchingNode,ProgressTimer>  itBranching (graph.iterator<BranchingNode>(), "statistics");

    // We loop the branching nodes
    for (itBranching.first(); !itBranching.isDone(); itBranching.next())
    {
        // We skip already visited nodes.
        if (marker.isMarked (*itBranching))  { continue; }

        // We launch the breadth first search; we get as a result the number of branching nodes in this component
        u_int64_t res = bfs.run (*itBranching);

        // We update our distribution
        distrib[res] ++;
    }

    // We compute the total number of branching nodes in all connected components.
    size_t sum = 0;   for (auto it = distrib.begin(); it != distrib.end(); it++)  {  sum += it->first*it->second; }

    // Note: it must be equal to the number of branching nodes of the graph
    assert (sum == itBranching.size());

    // We aggregate the computed information
    Properties props ("connected_components");
    props.add (1, "graph_name",              "%s", graph.getName().c_str());
    props.add (1, "nb_branching_nodes",      "%d", sum);
    props.add (1, "nb_connected_components", "%d", distrib.size());
    for (map<size_t,size_t>::iterator it = distrib.begin(); it!=distrib.end(); it++)
    {
        props.add (2, "component");
        props.add (3, "nb_nodes",    "%d", it->first);
        props.add (3, "nb_occurs",   "%d", it->second);
        props.add (3, "freq_nodes",  "%f", 100.0*(float)(it->first*it->second) / (float)sum);
        props.add (3, "freq_occurs", "%f", 100.0*(float)it->second / (float)sum);
    }

    // We dump the results in a XML file in the current directory
    XmlDumpPropertiesVisitor v (graph.getName() + ".xml", false);
    props.accept (&v);

    return EXIT_SUCCESS;
}
//! [snippet1]
