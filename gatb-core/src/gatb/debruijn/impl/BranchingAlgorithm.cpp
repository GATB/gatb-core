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

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace debruijn  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Graph: build branching nodes           ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
BranchingAlgorithm<span>::BranchingAlgorithm (
    const Graph& graph,
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
template<size_t span>
BranchingAlgorithm<span>::BranchingAlgorithm (tools::storage::impl::Storage& storage)
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
template<size_t span>
BranchingAlgorithm<span>::~BranchingAlgorithm ()
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
#ifndef WITH_LAMBDA_EXPRESSION
template<typename Count, typename Type>
struct FunctorNodes
{
    const Graph* graph; ThreadObject<u_int64_t>& count; ThreadObject <vector<Count> >& branchingNodes;
    FunctorNodes (const Graph* graph, ThreadObject<u_int64_t>& count, ThreadObject <vector<Count> >& branchingNodes)
        : graph(graph), count(count), branchingNodes(branchingNodes)  {}
    void operator() (const Node& node)
    {
        if (graph->isBranching(node))
        {
            count() ++;
            branchingNodes().push_back (Count (node.kmer.get<Type>(), node.abundance));
        }
    }
};
#endif
/*********************************************************************/

template<size_t span>
void BranchingAlgorithm<span>::execute ()
{
    /** We get a synchronized cache on the branching bag to be built. */
    ThreadObject<BagCache<Count> > branching = BagCache<Count> (_branchingCollection, 8*1024, System::thread().newSynchronizer());
    ThreadObject<u_int64_t> count;

    /** We get an iterator over all graph nodes. */
    Graph::Iterator<Node> itNodes = _graph->iterator<Node>();

    /** We encapsulate this iterator with a potentially decorated iterated (for progress information). */
    tools::dp::Iterator<Node>* iter = createIterator<Node> (
        itNodes.get(),
        itNodes.size(),
        progressFormat1
    );
    LOCAL (iter);

    ThreadObject <vector<Count> > branchingNodes;

#ifdef WITH_LAMBDA_EXPRESSION
    auto functorNodes = [&] (const Node& node)
    {
        if (this->_graph->isBranching(node))
        {
            count() ++;
            branchingNodes().push_back (Count (node.kmer.get<Type>(), node.abundance));
        }
    };
#else
    FunctorNodes<Count,Type> functorNodes (this->_graph, count, branchingNodes);
#endif

    /** We iterate the nodes. */
    tools::dp::IDispatcher::Status status = getDispatcher()->iterate (iter, functorNodes);

    /** We concate the kmers. */
    for (size_t i=0; i<branchingNodes.size(); i++)
    {
        vector<Count>& v = branchingNodes[i];
        branchingNodes->insert (branchingNodes->end(), v.begin(), v.end());
        v.clear ();
    };

    /** We sort the kmers. */
    sort (branchingNodes->begin(), branchingNodes->end());

    /** We put the kmers into the final bag. */
    _branchingCollection->insert (branchingNodes->data(), branchingNodes->size());

    /** We save the kind in the storage. */
    _storage(getName()).addProperty ("kind", toString(_kind));

    /** We gather the information collected during iteration. */
    for (size_t i=0; i<count.size(); i++)  {  *count += count[i];  }

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "nb_branching", "%ld", *count);
    getInfo()->add (2, "percentage",   "%.1f", (itNodes.size() > 0 ? 100.0*(float)*count/(float)itNodes.size() : 0));

    getInfo()->add (1, "time");
    getInfo()->add (2, "build", "%.3f", status.time / 1000.0);
}


/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class BranchingAlgorithm <KSIZE_1>;
template class BranchingAlgorithm <KSIZE_2>;
template class BranchingAlgorithm <KSIZE_3>;
template class BranchingAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
