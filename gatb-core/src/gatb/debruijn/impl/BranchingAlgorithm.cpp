/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
template<typename ProductFactory, size_t span>
BranchingAlgorithm<ProductFactory,span>::BranchingAlgorithm (
    const Graph& graph,
    tools::collections::Bag<Count>* branchingBag,
    size_t                      nb_cores,
    tools::misc::IProperties*   options
)
    : Algorithm("branching", nb_cores, options), _graph (graph),  _branchingBag(0)
{
    setBranchingBag(branchingBag);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, size_t span>
BranchingAlgorithm<ProductFactory,span>::~BranchingAlgorithm ()
{
    setBranchingBag(0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename ProductFactory, size_t span>
void BranchingAlgorithm<ProductFactory,span>::execute ()
{
    /** We get a synchronized cache on the branching bag to be built. */
    ThreadObject<BagCache<Count> > branching = BagCache<Count> (_branchingBag, 8*1024, System::thread().newSynchronizer());
    ThreadObject<u_int64_t> count;

    /** We get an iterator over all graph nodes. */
    Graph::Iterator<Node> itNodes = _graph.iterator<Node>();

    /** We encapsulate this iterator with a potentially decorated iterated (for progress information). */
    tools::dp::Iterator<Node>* iter = createIterator<Node> (
        itNodes.get(),
        itNodes.size(),
        progressFormat1
    );
    LOCAL (iter);

    ThreadObject <vector<Count> > branchingNodes;

    /** We iterate the nodes. */
    tools::dp::IDispatcher::Status status = getDispatcher()->iterate (iter, [&] (const Node& node)
    {
        if (this->_graph.isBranching(node))
        {
            count() ++;
            branchingNodes().push_back (Count (node.kmer.get<Type>(), node.abundance));
        }
    });

    /** We concate the kmers. */
    branchingNodes.foreach ([&] (vector<Count>& v)
    {
        branchingNodes->insert (branchingNodes->end(), v.begin(), v.end());
        v.clear ();
    });

    /** We sort the kmers. */
    sort (branchingNodes->begin(), branchingNodes->end());

    /** We put the kmers into the final bag. */
    _branchingBag->insert (branchingNodes->data(), branchingNodes->size());

    /** We gather the information collected during iteration. */
    count.foreach    ([&] (u_int64_t n) { *count += n; } );

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

template class BranchingAlgorithm <ProductFileFactory, 32*1>;
template class BranchingAlgorithm <ProductFileFactory, 32*2>;
template class BranchingAlgorithm <ProductFileFactory, 32*3>;
template class BranchingAlgorithm <ProductFileFactory, 32*4>;

/********************************************************************************/

template class BranchingAlgorithm <ProductHDF5Factory, 32*1>;
template class BranchingAlgorithm <ProductHDF5Factory, 32*2>;
template class BranchingAlgorithm <ProductHDF5Factory, 32*3>;
template class BranchingAlgorithm <ProductHDF5Factory, 32*4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
