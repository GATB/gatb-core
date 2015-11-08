
#ifndef _GATB_CORE_DEBRUIJN_IMPL_INSTANTIATIONS_
#define _GATB_CORE_DEBRUIJN_IMPL_INSTANTIATIONS_


#include <gatb/tools/math/Integer.hpp>

#include <gatb/debruijn/impl/Simplifications.cpp> // important that it's .cpp and not .hpp, see https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
#include <gatb/debruijn/impl/Graph.cpp> 
#include <gatb/debruijn/impl/Terminator.cpp> 
#include <gatb/debruijn/impl/Traversal.cpp> 
#include <gatb/debruijn/impl/Frontline.cpp> 

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

typedef GraphTemplate<Node, Edge, GraphDataVariant> GraphT;
typedef Node NodeT;
typedef Edge EdgeT;
typedef BranchingNode_t<Node> BranchingNodeT;
typedef BranchingEdge_t<Node, Edge> BranchingEdgeT;

#include <gatb/debruijn/impl/Instantiations.hpp> // this might be a bit exotic to do it like that.. but it works.

// uses Node and Edge as defined in Graph.hpp (legacy GATB compatibility, when Graph was not templated)
template class GraphTemplate<Node, Edge, GraphDataVariant>; 

template class Simplifications<Node, Edge, GraphDataVariant>; 

// leagcy GATB compatibility
template class MPHFTerminatorTemplate<Node, Edge, GraphDataVariant>; 
template class BranchingTerminatorTemplate<Node, Edge, GraphDataVariant>; 
template class TerminatorTemplate<Node, Edge, GraphDataVariant>; 

// legacy GATB compatibility
template class TraversalTemplate<Node, Edge, GraphDataVariant>; 
template class MonumentTraversalTemplate<Node, Edge, GraphDataVariant>; 
template class SimplePathsTraversalTemplate<Node, Edge, GraphDataVariant>; 
template class NullTraversalTemplate<Node, Edge, GraphDataVariant>; 

// legacy GATB compatbility
template class FrontlineTemplate<Node, Edge, GraphDataVariant>; 
template class FrontlineBranchingTemplate<Node, Edge, GraphDataVariant>; 
template class FrontlineReachableTemplate<Node, Edge, GraphDataVariant>; 



/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
