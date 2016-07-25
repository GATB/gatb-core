/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Uricaru, P.Peterlongo, R.Chikhi, G.Rizk, E.Drezen
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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_ITERATIVE_EXTENSION_H
#define _GATB_CORE_DEBRUIJN_IMPL_ITERATIVE_EXTENSION_H

/********************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/Terminator.hpp>
#include <gatb/debruijn/impl/Traversal.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

/** \brief Class providing helpers for graph traversal
 *
 * The implementation relies on the Graph class of GATB-CORE, which provides the
 * API to traverse a de Bruijn graph.
 */
template <size_t span=KMER_DEFAULT_SPAN, typename Node=Node_t<>, typename Edge=Edge_t<Node_t<> >, typename Graph_t=Graph>
class IterativeExtensions
{
public:

    /** ShortcutS. */
    typedef typename kmer::impl::Kmer<span>::ModelCanonical       Model;
    typedef typename kmer::impl::Kmer<span>::ModelCanonical::Kmer KmerModel;
    typedef typename kmer::impl::Kmer<span>::Type                 kmer_type;

    /** Constructor.
     * \param[in] graph : de Bruijn graph of the reads
     * \param[in] terminator : object used to mark nodes during graph traversal.
     * \param[in] traversalKind : kind of graph traversal (unitig,contig)
     * \param[in] whenToStop : tells when to stop the extension
     * \param[in] searchMode : mode of graph traversal
     * \param[in] dontOutputFirstNucl : tells whether the first nucleotide of an extension has to be output
     * \param[in] max_depth : max depth of the graph traversal
     * \param[in] max_nodes : max nodes
     */
    IterativeExtensions (
        const Graph_t&                  graph,
        TerminatorTemplate<Node,Edge,Graph_t>&                   terminator,
        tools::misc::TraversalKind    traversalKind,
        tools::misc::ExtendStopMode_e whenToStop,
        tools::misc::SearchMode_e     searchMode,
        bool                          dontOutputFirstNucl,
        int                           max_depth,
        int                           max_nodes
    );

    /** Return the unitig/contig which starts with L, as well as all ths contigs that follow him, up to max_depth.
     * Results go to a fasta file (output_file)
     * \param[in] L : starter string for the unitig/contig
     * \param[in] R :
     * \param[in] outputBank : bank where the unitig/contig have to be dumped
     * \param[in] swf : to be detailed
     */
    void construct_linear_seqs (
        const std::string& L,
        const std::string& R,
        bank::IBank*       outputBank,
        bool               swf = false
    );

    /** Return the unitig/contig which starts with L, as well as all ths contigs that follow him, up to max_depth.
     * Results go to a fasta file (output_file)
     * \param[in] L : starter string for the unitig/contig
     * \param[in] R :
     * \param[in] output_file : uri of the file where the unitig/contig have to be dumped
     * \param[in] swf : to be detailed
     */
    void construct_linear_seqs (
        const std::string& L,
        const std::string& R,
        const std::string& output_file,
        bool               swf = false
    );

private:

    const Graph_t&                  graph;
    TerminatorTemplate<Node,Edge,Graph_t>&                   terminator;
    tools::misc::TraversalKind      traversalKind;
    tools::misc::ExtendStopMode_e   when_to_stop_extending;
    tools::misc::SearchMode_e       searchMode;
    bool                            dont_output_first_nucleotide;
    int                             max_depth;
    int                             max_nodes;

    Model model;
    Model modelMinusOne;

    /** Fill a Sequence instance from the results of the current graph traversal.
     * \param[in] node : starting node of the path
     * \param[in] consensusRight : consensus nucleotides to be assign to the Sequence instance
     * \param[in] nbNodes : nb nodes
     * \param[in] depth : depth
     * \param[out] seq : the sequence to be filled.
     */
    void buildSequence (
        const Node&     node,
        const Path_t<Node>&     consensusRight,
        size_t          nbNodes,
        size_t          depth,
        bank::Sequence& seq
    );

    /** */
    bool compare_and_mark_last_k_minus_one_mer (const std::string& node, std::set<kmer_type>& kmers_set);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_ITERATIVE_EXTENSION_H */
