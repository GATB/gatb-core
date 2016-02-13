//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <map>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

using namespace std;

/********************************************************************************/
/*                   Statistics about branching nodes in/out degrees.           */
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
    parser.push_back (new OptionOneParam (STR_NB_CORES,  "nb cores",     false, "0"));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        // We load the graph
        Graph graph = Graph::load (options->getStr(STR_URI_GRAPH));

        // We set the number of cores to be used. Use all available cores if set to 0.
        size_t nbCores = options->getInt(STR_NB_CORES);

        // We get an iterator for branching nodes of the graph.
        // We use a progress iterator to get some progress feedback
        ProgressGraphIterator<BranchingNode,ProgressTimer>  itBranching (graph.iteratorBranching(), "statistics");

        // We define some kind of unique identifier for a couple (indegree,outdegree)
        typedef pair<size_t,size_t> InOut_t;

        // We want to gather some statistics during the iteration.
        // Note the use of ThreadObject: this object will be cloned N times (one object per thread) and each clone will
        // be reachable within the iteration block through ThreadObject::operator()
        ThreadObject <map <InOut_t, size_t> > topology;

        // We dispatch the iteration on several cores. Note the usage of lambda expression here.
        IDispatcher::Status status = Dispatcher(nbCores).iterate (itBranching, [&] (const BranchingNode& node)
        {
            // We retrieve the current instance of map <InOut_t,size_t> for the current running thread.
            map <InOut_t,size_t>& localTopology = topology();

            // We get branching nodes neighbors for the current branching node.
            Graph::Vector<BranchingEdge> successors   = graph.successorsBranchingEdge ((Node&)node); /* node casting is necessary due to incomplete Graph.hpp api for successorsBranchingEdge */
            Graph::Vector<BranchingEdge> predecessors = graph.predecessorsBranchingEdge ((Node&)node); 

            // We increase the occurrences number for the current couple (in/out) neighbors
            localTopology [make_pair(predecessors.size(), successors.size())] ++;
        });

        // Now, the parallel processing is done. We want now to aggregate the information retrieved
        // in each thread in a single map.

        // We get each map<InOut_t,size_t> object filled in each thread, and we add its data into the "global" map.
        // The global map is reachable through the ThreadObject::operator*. The "topology.foreach" will loop over
        // all cloned object used in the threads.
        topology.foreach ([&] (const map <InOut_t, size_t>& t)
        {
            // We update the occurrence of the current couple (in/out)
            for_each (t.begin(), t.end(), [&] (const pair<InOut_t, size_t>& p) { (*topology)[p.first] += p.second;  });
        });

        // We sort the statistics by decreasing occurrence numbers. Since map have its own ordering, we need to put all
        // the data into a vector and sort it with our own sorting criteria.
        vector < pair<InOut_t,size_t> >  stats;
        for (auto it = topology->begin(); it != topology->end(); it++)  { stats.push_back (*it); }
        sort (stats.begin(), stats.end(), [=] (const pair<InOut_t,size_t>& a, const pair<InOut_t,size_t>& b) { return a.second > b.second; });

        printf ("\nThere are %d branching nodes with the following distribution: \n", itBranching.size());

        size_t sum=0;
        for (size_t i=0; i<stats.size(); i++)
        {
            sum += stats[i].second;

            printf ("    [in=%d out=%d]  nb=%7d  percent=%5.2f  distrib=%5.2f\n",
                stats[i].first.first,
                stats[i].first.second,
                stats[i].second,
                100.0*(float)stats[i].second / (float)itBranching.size(),
                100.0*(float)sum             / (float)itBranching.size()
            );
        }

        printf ("\nDone on %d cores in %.2f sec\n\n", status.nbCores, (float)status.time/1000.0);
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
