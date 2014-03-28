/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

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

// We define a class that marks nodes in a graph and call tell whether a given node is marked or not.
// We use a map for the implementation (could be not optimal).
template<typename T>  class GraphMarker
{
public:

    GraphMarker (const Graph& graph) : graph(graph)
    {
        // We insert all the nodes into our map.
        graph.iterator<T>().iterate (IterateNodes(this));
    }

    void mark (const T& item)  {  markMap [item] = true;  }

    void mark (const set<T>& items)
    {
        for (typename set<T>::const_iterator it=items.begin(); it != items.end(); ++it)  {  mark (*it); }
    }

    bool isMarked (const T& item) const  {  return markMap.find(item)->second;  }

private:
    const Graph& graph;
    map<T,bool>  markMap;

    struct IterateNodes
    {
        IterateNodes (GraphMarker* marker) : marker(marker) {}
        GraphMarker* marker;
        void operator() (const T& item) const  {  marker->markMap[item] = false;  }
    };
};

/********************************************************************************/

// We define a class that performs a breadth first search in a graph from a starting node.
// We need a marker to globally mark the visited nodes from one call to another.
// The nodes visit is handled by a std::queue object.
template<typename T>  class BFS
{
public:

    BFS (const Graph& graph) : graph(graph) {}

    // Process the BFS started from the provided node. As a result, we get the nodes number in the
    // found connected component.
    const set<T>& run (const T& node)
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
        set<T> frontline;
        frontline.insert (node);

        // We set the initial state for C(0)
        connectedComponent.clear();
        connectedComponent.insert (node);

        // We clear the statistics
        bfsInfo.clear ();

        // We launch the recursion.
        while (!frontline.empty())
        {
            // We add statistics information for the current depth in the BFS.
            bfsInfo.push_back (make_pair(frontline.size(), connectedComponent.size()));

            // We get the neighbors for the current front line, ie. we get N(F(n))
            set<T> neighbors = graph.neighbors<T> (frontline.begin(), frontline.end());

            // We reset the current front line => we reuse it for computing F(n+1)
            frontline.clear();

            // We compute the recursion for F(n+1) and C(n+1)
            for (typename set<T>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
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
    const set<T>& get() const { return connectedComponent; }

    // Statistics
    const list<pair<size_t,size_t> >& getStatistics() const { return bfsInfo; }

private:
    const Graph& graph;
    set<T>       connectedComponent;

    list<pair<size_t,size_t> > bfsInfo;
};

/** */
struct Entry   {  size_t nbOccurs;  size_t nbKmers;    Entry() : nbOccurs(0), nbKmers(0) {} };

/** */
typedef pair<size_t,size_t> InOut_t;
bool CompareFct (const pair<InOut_t,size_t>& a, const pair<InOut_t,size_t>& b) { return a.second > b.second; }

/********************************************************************************/
/*        Computing connected components of the branching nodes subgraph.       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.push_back (new OptionOneParam (STR_URI_INPUT,  "graph file", true));

    IProperties* params = 0;

    try  {
        /** We parse the user options. */
        params = parser.parse (argc, argv);
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    // We create the graph with the bank and other options
    Graph graph = Graph::load (params->getStr(STR_URI_INPUT));

    // We create a graph marker.
    GraphMarker<BranchingNode> marker (graph);

    // We create an object for Breadth First Search for the de Bruijn graph.
    BFS<BranchingNode> bfs (graph);

    // We want to compute the distribution of connected components of the branching nodes.
    //    - key is a connected component class (for a given number of branching nodes for this component)
    //    - value is the number of times this component class occurs in the branching sub graph
    map<size_t,Entry> distrib;

    // We get an iterator for all nodes of the graph. We use a progress iterator to get some progress feedback
    ProgressGraphIterator<BranchingNode,ProgressTimer>  itBranching (graph.iterator<BranchingNode>(), "statistics");

    // We want to know the number of connected components
    size_t nbConnectedComponents = 0;

    // We define some kind of unique identifier for a couple (indegree,outdegree)
    map <InOut_t, size_t> topology;

    size_t simplePathSizeMin = ~0;
    size_t simplePathSizeMax =  0;


    // We want time duration of the iteration
    TimeInfo ti;
    ti.start ("compute");

    // We loop the branching nodes
    for (itBranching.first(); !itBranching.isDone(); itBranching.next())
    {
        // We get branching nodes neighbors for the current branching node.
        Graph::Vector<BranchingEdge> successors   = graph.successors  <BranchingEdge> (*itBranching);
        Graph::Vector<BranchingEdge> predecessors = graph.predecessors<BranchingEdge> (*itBranching);

        // We increase the occurrences number for the current couple (in/out) neighbors
        topology [make_pair(predecessors.size(), successors.size())] ++;

        // We loop the in/out neighbors and update min/max simple path size
        for (size_t i=0; i<successors.size(); i++)
        {
            simplePathSizeMax = std::max (simplePathSizeMax, successors[i].distance);
            simplePathSizeMin = std::min (simplePathSizeMin, successors[i].distance);
        }
        for (size_t i=0; i<predecessors.size(); i++)
        {
            simplePathSizeMax = std::max (simplePathSizeMax, predecessors[i].distance);
            simplePathSizeMin = std::min (simplePathSizeMin, predecessors[i].distance);
        }

        // We skip already visited nodes.
        if (marker.isMarked (*itBranching))  { continue; }

        // We launch the breadth first search; we get as a result the set of branching nodes in this component
        const set<BranchingNode>& component = bfs.run (*itBranching);

        // We mark the nodes for this connected component
        marker.mark (component);

        // We update our distribution
        distrib[component.size()].nbOccurs += 1;

        // We update the number of connected components.
        nbConnectedComponents++;
    }

    ti.stop ("compute");

    // We compute the total number of branching nodes in all connected components.
    size_t sumOccurs = 0;
    size_t sumKmers = 0;
    for (map<size_t,Entry>::iterator it = distrib.begin(); it != distrib.end(); it++)
    {
        sumOccurs += it->first*it->second.nbOccurs;
        sumKmers  += it->second.nbKmers;
    }

    // We sort the statistics by decreasing occurrence numbers. Since map have its own ordering, we need to put all
    // the data into a vector and sort it with our own sorting criteria.
    vector < pair<InOut_t,size_t> >  stats;
    for (map <InOut_t, size_t>::iterator it = topology.begin(); it != topology.end(); it++)  { stats.push_back (*it); }

    sort (stats.begin(), stats.end(), CompareFct);

    // Note: it must be equal to the number of branching nodes of the graph
    assert (sumOccurs == itBranching.size());

    // We aggregate the computed information
    Properties props ("topology");

    props.add (1, "graph");
    props.add (2, "name",                    "%s", graph.getName().c_str());
    props.add (2, "db_input",                "%s", graph.getInfo().getStr("input").c_str());
    props.add (2, "db_nb_seq",               "%d", graph.getInfo().getInt("sequences_number"));
    props.add (2, "db_size",                 "%d", graph.getInfo().getInt("sequences_size"));
    props.add (2, "kmer_size",               "%d", graph.getInfo().getInt("kmer_size"));
    props.add (2, "kmer_nks",                "%d", graph.getInfo().getInt("nks"));
    props.add (2, "nb_nodes",                "%d", graph.getInfo().getInt("kmers_nb_solid"));
    props.add (2, "nb_branching_nodes",      "%d", graph.getInfo().getInt("nb_branching"));
    props.add (2, "percent_branching_nodes", "%.1f",
        graph.getInfo().getInt("kmers_nb_solid") > 0 ?
        100.0 * (float)graph.getInfo().getInt("nb_branching") / (float) graph.getInfo().getInt("kmers_nb_solid") : 0
    );

    props.add (1, "branching_nodes");

    props.add (2, "simple_path");
    props.add (3, "size_min", "%d", simplePathSizeMin);
    props.add (3, "size_max", "%d", simplePathSizeMax);

    props.add (2, "neighborhoods");
    for (size_t i=0; i<stats.size(); i++)
    {
        props.add (3, "neighborhood", "in=%d out=%d", stats[i].first.first, stats[i].first.second);
        props.add (4, "nb_bnodes",     "%d",    stats[i].second);
        props.add (4, "percentage",   "%5.2f", itBranching.size() > 0 ?
            100.0*(float)stats[i].second / (float)itBranching.size() : 0
        );
    }

    props.add (2, "connected_components");
    props.add (3, "nb_classes",    "%d", distrib.size());
    props.add (3, "nb_components", "%d", nbConnectedComponents);
    for (map<size_t,Entry>::iterator it = distrib.begin(); it!=distrib.end(); it++)
    {
        props.add (3, "component_class");
        props.add (4, "nb_occurs",    "%d", it->second.nbOccurs);
        props.add (4, "nb_bnodes",    "%d", it->first);
        props.add (4, "freq_bnodes",  "%f", sumOccurs > 0 ?
            100.0*(float)(it->first*it->second.nbOccurs) / (float)sumOccurs : 0
        );
    }
    props.add (1, ti.getProperties("time"));

    // We dump the results in a XML file in the current directory
    XmlDumpPropertiesVisitor v (graph.getName() + ".xml", false);
    props.accept (&v);

    return EXIT_SUCCESS;
}
//! [snippet1]
