//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <set>

#undef NDEBUG
#include <cassert>

using namespace std;

/********************************************************************************/

//
//
/********************************************************************************/
/*  BFS exploration of a graph.                                                 */
/*                                                                              */
/*  We define a class that marks nodes in a graph and call tell whether a given */
/*  node is marked or not.                                                      */
/*  We use a map for the implementation (could be not optimal).                 */
/*                                                                              */
/* Cmd-line: debruijn19 -graph <h5 file>                                        */
/*                                                                              */
/* Sample: debruijn19 -graph gatb-core/gatb-core/test/db/celegans_reads.h5      */
/*                                                                              */
/* Note:                                                                        */
/*     - '.h5' file contains the HDF5 formatted representation of a de bruijn   */
/*     graph created from a set of reads.                                       */
/*     - a '.h5' file is created using dbgh5 program provided with GATB-Core.   */
/*     Basic use is as follows:                                                 */
/*        dbgh5 -in <fasta/q file> -out <h5 file>                               */
/*     You can also control kmer-size and kmer abundance, see dbgh5 help.       */
/*                                                                              */
/********************************************************************************/
class GraphMarker
{
public:

    GraphMarker (const Graph& graph) : graph(graph)
    {
        // We insert all the nodes into our map.
        auto iter = [&] (const BranchingNode& item)  {  this->markMap[item] = false;  };
        graph.iteratorBranching().iterate (iter);
    }

    void mark (const BranchingNode& item)  {  markMap [item] = true;  }

    void mark (const set<BranchingNode>& items)
    {
        for (typename set<BranchingNode>::const_iterator it=items.begin(); it != items.end(); ++it)  {  mark (*it); }
    }

    bool isMarked (const BranchingNode& item) const  {  return markMap.find(item)->second;  }

private:
    const Graph& graph;
    map<BranchingNode,bool>  markMap;
};

/********************************************************************************/

// We define a class that performs a breadth first search in a graph from a starting node.
// We need a marker to globally mark the visited nodes from one call to another.
// The nodes visit is handled by a std::queue object.
class BFS
{
public:

    BFS (const Graph& graph) : graph(graph) {}

    // Process the BFS started from the provided node. As a result, we get the nodes number in the
    // found connected component.
    const set<BranchingNode>& run (const BranchingNode& node)
    {
        // ALGORITHM (recursion on n)
        //    We define as C(n) the set of nodes of the connected component to be computed.
        //    We define as F(n) the set of nodes for which we get neighbors for extending the BFS
        //    We use the function N(E) that returns the set of neighbors nodes of items in set E
        //
        //  Init:
        //      C(0)={b} and F(0)={b},  where b is the initial node
        //  Recursion:
        //      F(n+1) = N(F(n)) - C(n)
        //      C(n+1) = N(F(n)) + C(n)
        //  End criteria:
        //      card(F(n)) = 0

        // We set the initial state for F(0)
        set<BranchingNode> frontline;
        frontline.insert (node);

        // We set the initial state for C(0)
        connectedComponent.clear();
        connectedComponent.insert (node);

        // We launch the recursion.
        while (!frontline.empty())
        {
            // We get the neighbors for the current front line, ie. we get N(F(n))
            set<BranchingNode> neighbors = graph.neighbors (frontline.begin(), frontline.end()); // this is a very special overload of neighbors() from Graph.hpp

            // We reset the current front line => we reuse it for computing F(n+1)
            frontline.clear();

            // We compute the recursion for F(n+1) and C(n+1)
            for (typename set<BranchingNode>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
            {
                if (connectedComponent.find (*it) == connectedComponent.end())
                {
                    // F(n+1) = N(F(n)) - C(n)
                    frontline.insert (*it);

                    // C(n+1) = N(F(n)) + C(n)
                    connectedComponent.insert (*it);
                }
            }
        }

        // We return the number of nodes for this connected component
        return connectedComponent;
    }

    // We provide an accessor to the nodes of the found connected component
    const set<BranchingNode>& get() const { return connectedComponent; }

private:
    const Graph& graph;
    set<BranchingNode>       connectedComponent;
};

/********************************************************************************/
/*        Computing connected components of the branching nodes subgraph.       */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("GraphStats");
    parser.push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",  true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We load the graph
        Graph graph = Graph::load (options->getStr(STR_URI_GRAPH));

        // We create a graph marker.
        GraphMarker marker (graph);

        // We create an object for Breadth First Search for the de Bruijn graph.
        BFS bfs (graph);

        // We want to compute the distribution of connected components of the branching nodes.
        //    - key is a connected component class (for a given number of branching nodes for this component)
        //    - value is the number of times this component class occurs in the branching sub graph
        map<size_t,size_t> distrib;

        // We get an iterator for branching nodes of the graph. We use a progress iterator to get some progress feedback
        ProgressGraphIterator<BranchingNode,ProgressTimer>  itBranching (graph.iteratorBranching(), "statistics");

        // We want to know the number of connected components
        size_t nbConnectedComponents = 0;

        // We want time duration of the iteration
        TimeInfo ti;
        ti.start ("compute");

        // We loop the branching nodes
        for (itBranching.first(); !itBranching.isDone(); itBranching.next())
        {
            // We skip already visited nodes.
            if (marker.isMarked (*itBranching))  { continue; }

            // We launch the breadth first search; we get as a result the set of branching nodes in this component
            const set<BranchingNode>& component = bfs.run (*itBranching);

            // We mark the nodes for this connected component
            marker.mark (component);

            // We update our distribution
            distrib[component.size()] ++;

            // We update the number of connected components.
            nbConnectedComponents++;
        }

        ti.stop ("compute");

        // We compute the total number of branching nodes in all connected components.
        size_t sum = 0;   for (map<size_t,size_t>::iterator it = distrib.begin(); it != distrib.end(); it++)  {  sum += it->first*it->second; }

        // Note: it must be equal to the number of branching nodes of the graph
        assert (sum == itBranching.size());

        // We aggregate the computed information
        Properties props ("connected_components");
        props.add (1, "graph_name",              "%s", graph.getName().c_str());
        props.add (1, "nb_branching_nodes",      "%d", sum);
        props.add (1, "nb_connected_classes",    "%d", distrib.size());
        props.add (1, "nb_connected_components", "%d", nbConnectedComponents);
        for (map<size_t,size_t>::iterator it = distrib.begin(); it!=distrib.end(); it++)
        {
            props.add (2, "component");
            props.add (3, "nb_nodes",    "%d", it->first);
            props.add (3, "nb_occurs",   "%d", it->second);
            props.add (3, "freq_nodes",  "%f", 100.0*(float)(it->first*it->second) / (float)sum);
            //props.add (3, "freq_occurs", "%f", 100.0*(float)it->second / (float)sum);
        }
        props.add (1, ti.getProperties("time"));

        // We dump the results in a XML file in the current directory
        XmlDumpPropertiesVisitor v (graph.getName() + ".xml", false);
        props.accept (&v);
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
