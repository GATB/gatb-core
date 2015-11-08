
#ifndef _GATB_CORE_DEBRUIJN_IMPL_INSTANTIATIONS_
#define _GATB_CORE_DEBRUIJN_IMPL_INSTANTIATIONS_


#include <gatb/tools/math/Integer.hpp>

#include <gatb/debruijn/impl/Graph.hpp>

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

typedef GraphTemplate<Node, Edge, GraphDataVariant> GraphT;
typedef Node NodeT;
typedef Edge EdgeT;
typedef BranchingNode_t<Node> BranchingNodeT;
typedef BranchingEdge_t<Node, Edge> BranchingEdgeT;

#include <gatb/debruijn/impl/Instantiations.hpp> // this might be a bit exotic to do it like that.. but it works.


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
