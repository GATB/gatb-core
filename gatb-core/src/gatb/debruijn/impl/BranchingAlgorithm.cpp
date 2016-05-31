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

#include <gatb/debruijn/impl/BranchingAlgorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <queue>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace debruijn  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Graph: build branching nodes           ";
static const char* progressFormat2 = "Graph: nb branching found : %-9d  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Node, typename Edge, typename GraphDataVariant>
BranchingAlgorithm<span, Node, Edge, GraphDataVariant>::BranchingAlgorithm (
    const GraphTemplate<Node, Edge, GraphDataVariant>& graph,
    tools::storage::impl::Storage& storage,
    tools::misc::BranchingKind  kind,
    size_t                      nb_cores,
    tools::misc::IProperties*   options
)
    : Algorithm("branching", nb_cores, options), _graph (&graph),  _storage(storage), _kind(kind), _branchingCollection(0)
{
    setBranchingCollection (& storage("branching").getCollection<Count> ("nodes"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Node, typename Edge, typename GraphDataVariant>
BranchingAlgorithm<span, Node, Edge, GraphDataVariant>::BranchingAlgorithm (tools::storage::impl::Storage& storage)
    : Algorithm("branching", 0, 0), _graph(0), _storage(storage), _branchingCollection(0)
{
    setBranchingCollection (& storage("branching").getCollection<Count> ("nodes"));

    string xmlString = storage(this->getName()).getProperty ("xml");
    stringstream ss; ss << xmlString;   getInfo()->readXML (ss);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Node, typename Edge, typename GraphDataVariant>
BranchingAlgorithm<span, Node, Edge, GraphDataVariant>::~BranchingAlgorithm ()
{
    setBranchingCollection(0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Node, typename Edge, typename GraphDataVariant>
IOptionsParser* BranchingAlgorithm<span, Node, Edge, GraphDataVariant>::getOptionsParser ()
{
    IOptionsParser* parser = new OptionsParser ("branching");

    parser->push_back (new OptionOneParam (STR_BRANCHING_TYPE,    "branching type ('none' or 'stored')",      false, "stored"));
    parser->push_back (new OptionOneParam (STR_TOPOLOGY_STATS,    "topological information level (0 for none)", false, "0"));

    return parser;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
typedef pair<size_t,size_t> InOut_t;
inline bool CompareFct (const pair<InOut_t,size_t>& a, const pair<InOut_t,size_t>& b) { return a.second > b.second; }

template<typename Count, typename Type>
struct FunctorData
{
	FunctorData() {  }
    std::vector<Count> branchingNodes;
    map <InOut_t, size_t> topology;
};

template<typename Count, typename Type, typename Node, typename Edge, typename GraphDataVariant>
struct FunctorNodes
{
    const GraphTemplate<Node, Edge, GraphDataVariant>* graph;
    ThreadObject<FunctorData<Count,Type> >& functorData;

    FunctorNodes (const GraphTemplate<Node, Edge, GraphDataVariant>* graph, ThreadObject<FunctorData<Count,Type> >& functorData)
        : graph(graph), functorData(functorData)  {}

    void operator() (Node& node)
    {
        // We get branching nodes neighbors for the current node.
        typename GraphTemplate<Node, Edge, GraphDataVariant>::template Vector<Node> successors   = graph->template successors   (node);
        typename GraphTemplate<Node, Edge, GraphDataVariant>::template Vector<Node> predecessors = graph->template predecessors (node);

        if ( ! (successors.size()==1 && predecessors.size()==1) )
        {
            // the node is branching

        	FunctorData<Count,Type>& data = functorData();

        	data.branchingNodes.push_back (Count (node.template getKmer<Type>(), node.abundance));

        	data.topology [make_pair(predecessors.size(), successors.size())] ++;
        }
    }
};

/*********************************************************************/

template<typename Count>
class SortCmd : public tools::dp::ICommand, public system::SmartPointer
{
public:
    vector<Count>& vec;
    SortCmd (vector<Count>& vec) : vec(vec) {}
    void execute ()  {  std::sort (vec.begin(), vec.end());  }
};

/*********************************************************************/

template <typename T>
struct Compare
{
    /** Biggest values first. */
    bool operator() (const T& a, const T& b)  {  return ! (*(a.first) < *(b.first));  }
};

/*********************************************************************/

template<typename Count>
class CustomListener : public ProgressProxy
{
    Collection<Count>* _branchingCollection;

public:
    CustomListener (IteratorListener* ref, Collection<Count>* branchingCollection)
        : ProgressProxy(ref), _branchingCollection(branchingCollection) {}

    void finish ()  {}

    void finishPostponed ()
    {
        this->setMessage (Stringify::format(progressFormat2, _branchingCollection->getNbItems()));
        getRef()->finish();
    }
};

/*********************************************************************/

template <size_t span, typename Node, typename Edge, typename GraphDataVariant>
void BranchingAlgorithm<span, Node, Edge, GraphDataVariant>::execute ()
{
    /** We get an iterator over all graph nodes. */
    typename GraphTemplate<Node, Edge, GraphDataVariant>::template Iterator<Node> itNodes = _graph->GraphTemplate<Node, Edge, GraphDataVariant>::iterator();

    /** We create a custom listener that makes the finish() method, normally called at end of iteration, do nothing this time.
     * => we define our own 'finishPostponed' method that is called when all the information is ok. */
    CustomListener<Count>* listener = new CustomListener<Count> (
        createIteratorListener (itNodes.size(), progressFormat1),
        _branchingCollection
    );

    /** We encapsulate this iterator with a potentially decorated iterated (for progress information). */
    tools::dp::Iterator<Node>* iter = createIterator<Node> (
        itNodes.get(),
        itNodes.size(),
        progressFormat1,
        listener
    );
    LOCAL (iter);

    /** We get a synchronized object on the data handled by functors. */
    ThreadObject <FunctorData<Count,Type> > functorData;

    FunctorNodes<Count,Type, Node, Edge, GraphDataVariant> functorNodes (this->_graph, functorData);

    /** We iterate the nodes. */
    tools::dp::IDispatcher::Status status = getDispatcher()->iterate (iter, functorNodes);

    /** Now, because we iterated with N threads, we have N vector of branching nodes. (N=nbcores used by the dispatcher)
     *  We need to merge them.
     *  TODO doc: why do they need to be sorted actually?!
     *  The next step are:
     *      1) sort each vector
     *      2) sort/merge the N vectors in the final collection
     */

    /** Step 1 : sorting the N branching nodes vectors (with the dispatcher). */
    vector<tools::dp::ICommand*> sortCmds;
    for (size_t i=0; i<functorData.size(); i++)  {  sortCmds.push_back (new SortCmd<Count> (functorData[i].branchingNodes));  }
    getDispatcher()->dispatchCommands (sortCmds);

    /** Step 2 : sort/merge the N vectors.
     * We use a priority queue for the merge process. */
    typedef typename vector<Count>::iterator BranchingIterator;
    typedef pair<BranchingIterator,BranchingIterator> BranchingIteratorPair;

    /** We use a cache to improve IO performances. */
    CollectionCache<Count> branchingCache (*_branchingCollection, 16*1024, 0);

    Type checksum; checksum.setVal( 0);

    /** We initialize our priority queue. */
    priority_queue <BranchingIteratorPair, vector<BranchingIteratorPair>, Compare<BranchingIteratorPair> > pq;
    for (size_t i=0; i<functorData.size(); i++)
    {
        if (functorData[i].branchingNodes.empty() == false)
        {
            pq.push (make_pair (functorData[i].branchingNodes.begin(), functorData[i].branchingNodes.end()));
        }
    }

    /** We process the merge/sort. */
    while (!pq.empty())
    {
        /** We get the top iterators pair from the priority queue and remove from it. */
        BranchingIteratorPair it = pq.top();
        pq.pop();

        /** We check that the current iterator is not finished. */
        if (it.first != it.second)
        {
            /** We insert the Count object into the final bag. */
            branchingCache.insert (*it.first);

            /** Stats */
            checksum += it.first->value;

            /** We update the priority queue. */
            ++(it.first); if (it.first != it.second)  {  pq.push (it); }
        }
    }

    /** We have to flush the cache to be sure every items is put into the cached collection. */
    branchingCache.flush ();

    /** We call our 'custom' finish method. */
    listener->finishPostponed();

    /** We save the kind in the storage. */
    _storage(getName()).addProperty ("kind", toString(_kind));

    /* print the number of branching nodes (could be important for debugging, if a user experiences a crash and copypastes stdout) */
    //cout << "Graph has " << _branchingCollection->getNbItems() << " branching nodes." << endl;

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "nb_branching", "%ld", _branchingCollection->getNbItems());
    getInfo()->add (2, "percentage",   "%.1f", (itNodes.size() > 0 ? 100.0*(float)_branchingCollection->getNbItems()/(float)itNodes.size() : 0));

    stringstream ss;  ss << checksum;
    getInfo()->add (2, "checksum_branching", "%s", ss.str().c_str());

    if (getInput()->get(STR_TOPOLOGY_STATS) && getInput()->getInt(STR_TOPOLOGY_STATS) > 0)
    {
        /** We get some topological information. */
        for (size_t i=0; i<functorData.size(); i++)
        {
            for (map<InOut_t, size_t>::iterator it = functorData[i].topology.begin();  it != functorData[i].topology.end(); ++it)
            {
                functorData->topology[it->first] += it->second;
            }
        }

        /** We sort the statistics. */
        vector < pair<InOut_t,size_t> >  topologyStats;
        for (map<InOut_t, size_t>::iterator it = functorData->topology.begin();  it != functorData->topology.end(); ++it)  { topologyStats.push_back (*it); }
        sort (topologyStats.begin(), topologyStats.end(), CompareFct);

        getInfo()->add (1, "topology");
        for (size_t i=0; i<topologyStats.size(); i++)
        {
            getInfo()->add (3, "neighborhood", "[in=%ld out=%ld]  nb : %8ld   percent. : %5.2f",
                topologyStats[i].first.first, topologyStats[i].first.second, topologyStats[i].second,
                _branchingCollection->getNbItems() > 0 ?
                100.0*(float)topologyStats[i].second / (float)_branchingCollection->getNbItems() : 0
            );
        }
    }

    getInfo()->add (1, "time");
    getInfo()->add (2, "build", "%.3f", status.time / 1000.0);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
