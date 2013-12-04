/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Debloom, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file BranchingAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Algorithm that computes the branching nodes of a De Bruijn graph
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/tools/collections/impl/Product.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

template <typename ProductFactory, typename T>
class BranchingAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Constructor. */
    BranchingAlgorithm (
        const Graph& graph,
        tools::collections::Bag<kmer::Kmer<T> >* branchingBag,
        size_t                      nb_cores = 0,
        tools::misc::IProperties*   options  = 0
    );

    /** Destructor. */
    ~BranchingAlgorithm ();

    /** */
    void execute ();

private:

    const Graph& _graph;

    tools::collections::Bag<kmer::Kmer<T> >* _branchingBag;
    void setBranchingBag (tools::collections::Bag<kmer::Kmer<T> >* branchingBag)  {  SP_SETATTR(branchingBag); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_ */

