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
#include <gatb/tools/storage/impl/Storage.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

template <size_t span=KMER_DEFAULT_SPAN>
class BranchingAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::Model Model;
    typedef typename kmer::impl::Kmer<span>::Type  Type;
    typedef typename kmer::impl::Kmer<span>::Count Count;

    /** Constructor. */
    BranchingAlgorithm (
        const Graph& graph,
        tools::collections::Bag<Count>* branchingBag,
        size_t                      nb_cores = 0,
        tools::misc::IProperties*   options  = 0
    );

    /** Destructor. */
    ~BranchingAlgorithm ();

    /** */
    void execute ();

private:

    const Graph& _graph;

    tools::collections::Bag<Count>* _branchingBag;
    void setBranchingBag (tools::collections::Bag<Count>* branchingBag)  {  SP_SETATTR(branchingBag); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_BRANCHING_ALGORITHM_HPP_ */

