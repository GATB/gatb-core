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

#include <gatb/debruijn/impl/IterativeExtensions.hpp>
#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

#include <cstdio>
#define DEBUG(a)    //a
#define VERBOSE(a)
#define INFO(a)

using namespace std;

using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

//#define DONTMARK

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

// We define a structure holding a Node and a depth
template <typename Node>
struct NodeDepth
{
    Node node;
    int  depth;

    NodeDepth () : node(0,STRAND_FORWARD), depth(0) {}

    NodeDepth (typename Node::Value kmer, Strand strand, int depth) : node(kmer,strand), depth(depth) {}

    // needed for comparisons inside a list
    bool operator<(const NodeDepth &other) const
    {
        if (node.kmer  != other.node.kmer)  { return (node.kmer < other.node.kmer); }
        if (depth      != other.depth)      { return (depth     < other.depth);     }
        return (node.strand < other.node.strand);
    }
};

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
/*
 * our assembly graph is connected by (k-1)-overlaps,
 * so this function is used to make sure we see each (k-1)-overlap in at most one right extremity
 */
template<size_t span, typename Node, typename Edge, typename Graph>
bool IterativeExtensions<span, Node, Edge, Graph>::compare_and_mark_last_k_minus_one_mer (const string& node, set<kmer_type>& kmers_set)
{
    KmerModel leftKmer = modelMinusOne.codeSeed (node.c_str(), Data::ASCII, node.size() - modelMinusOne.getKmerSize());
    kmer_type kmer = leftKmer.value();

    if (kmers_set.find(kmer) != kmers_set.end())
        return true;

    kmers_set.insert(kmer);
    return false;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span, typename Node, typename Edge, typename Graph>
IterativeExtensions<span, Node, Edge, Graph>::IterativeExtensions (
    const Graph&        graph,
    TerminatorTemplate<Node,Edge,Graph>&         terminator,
    TraversalKind       traversalKind,
    ExtendStopMode_e    whenToStop,
    SearchMode_e        searchMode,
    bool                dontOutputFirstNucl,
    int                 max_depth,
    int                 max_nodes
)
    : graph(graph), terminator(terminator),
      traversalKind(traversalKind),  when_to_stop_extending(whenToStop),
      searchMode(searchMode), dont_output_first_nucleotide(dontOutputFirstNucl),
      max_depth(max_depth), max_nodes(max_nodes)
{
    model         = Model (graph.getKmerSize()  );
    modelMinusOne = Model (graph.getKmerSize()-1);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span, typename Node, typename Edge, typename Graph>
void IterativeExtensions<span, Node, Edge, Graph>::construct_linear_seqs (
    const string& L,
    const string& R,
    IBank*        outputBank,
    bool          swf
)
{
    /** Shortcuts. */
    size_t sizeKmer = graph.getKmerSize();

    /** We first reset the terminator. */
    terminator.reset(); // distinct extensions may share kmers, however, a unique extension doesn't.

    DEBUG ((cout << "[IterativeExtensions::construct_linear_seqs]  output=" << outputBank->getId() << "  L=" << L
        << "  search=" << searchMode << " swf=" << swf << endl
    ));

    /** We get a token on the bank. */
    LOCAL (outputBank);

    /** We create a Traversal instance. */
    TraversalTemplate<Node,Edge,Graph>* traversal = TraversalTemplate<Node,Edge,Graph>::create (traversalKind, graph, terminator, max_depth, 500, 20);
    LOCAL (traversal);

    long long nbNodes = 0;
    long long totalnt = 0;

    /** We need a container that holds NodeDepth objects during the extension. */
    vector <NodeDepth<Node> > kmers_to_traverse;

    /** We get the first kmer of the L string. */
    KmerModel leftKmer = model.codeSeed (L.c_str(), Data::ASCII, 0);

    /** We put this kmer into the vector of kmers to be processed. */
    NodeDepth<Node> ksd (typename Node::Value(leftKmer.value()), leftKmer.which() ? STRAND_FORWARD : STRAND_REVCOMP, 0);
    kmers_to_traverse.push_back (ksd);

    DEBUG ((cout << "---> kmer=" <<  leftKmer.value() << " strand=" << (leftKmer.which() ? "FW" : "RC") << endl));

#ifndef DONTMARK
    set<kmer_type> already_extended_from;
    //   compare_and_mark_last_k_minus_one_mer(L, already_extended_from); // mark first kmer to never extend from it again, // L will be marked at first iteration below
#endif

    /** We will need a Path object and a Sequence object during the extension. */
    Path_t<Node> rightTraversal;
    Node endNode;
    Sequence seq (Data::ASCII);

    /**************************************************************/
    /**          MAIN LOOP ON THE REMAINING KMERS                 */
    /**************************************************************/
    while (kmers_to_traverse.size() > 0) // min_depth is max_gap_length here
    {
        VERBOSE (("IterativeExtensions::construct_linear_seqs  MAIN LOOP %ld\n", kmers_to_traverse.size()));

        if (searchMode == SearchMode_Depth)
        {
            ksd = kmers_to_traverse.back();
            kmers_to_traverse.pop_back();
        }
        else if (searchMode == SearchMode_Breadth)
        {
            ksd = kmers_to_traverse.front();
            kmers_to_traverse.erase (kmers_to_traverse.begin());
        }

        /** We compute the extension on the right. */
        int len_right = traversal->traverse (ksd.node, endNode, DIR_OUTCOMING, rightTraversal);

        DEBUG ((cout << "------> kmer=" << std::hex << ksd.node.kmer.get<kmer_type>() << std::dec
            << "  strand=" << toString(ksd.node.strand) << "  depth=" << ksd.depth
            << "  len_right=" << len_right << endl
        ));

        /** We build the sequence to be inserted in the output bank. */
        buildSequence (ksd.node, rightTraversal, nbNodes, ksd.depth, seq);

        /** We insert the sequence into the output bank. */
        outputBank->insert (seq);

        /** We update statistics. */
        int node_len = len_right + sizeKmer;
        nbNodes += 1;
        totalnt += node_len;

        // if we only want 1 extension, stop now
        if (when_to_stop_extending == ExtendStopMode_after_first_contig)
        {
            INFO (("Stopping because we want only 1 extension\n"));
            break;
        }

        if (swf)
        {

	    /*  old version            
	    char* found = strstr (seq.getDataBuffer(), R.c_str());
            if (found != NULL  &&  ksd.depth > (int)sizeKmer)
            {
                INFO (("swf STOP \n"));
                break;
            }
 */
        	// Apr 2017 : new version for MindTheGap fill with contigs
        	//target can be the concatenation of several target kmers, checks for all the target kmers, if one found stop extending this contig but continue other branches
        	std::string subseed;
        	std::string target= R;
        	bool stopExtend= false;
        	for (unsigned i = 0; i < target.length(); i += sizeKmer)
        	{
        		subseed=target.substr(i, sizeKmer);
        		const char* found = strstr (seq.toString().c_str(), subseed.c_str());

        		if (found != NULL  &&  ksd.depth > (int)sizeKmer)
        		{
        			stopExtend=true;
        			break;
        		}

        	}
        	if (stopExtend) continue; //one of the targets was found, stop this extension but keep extending other branches

        }

        if (nbNodes > max_nodes) //GR stop when too complex  huum when to stop ?
        {
            INFO (("... XXX Stopped extending node %s because %lld nodes reached. Was at depth %d.\n",
                seq.toString().c_str(), nbNodes, ksd.depth
            ));
            break;
        }


        // if max depth reached, don't extend that one
        if (ksd.depth + node_len > max_depth)
        {
            INFO (("... XXX max depth reached for node %s (depth + node length %i %i = %i) \n",
                seq.toString().c_str(), ksd.depth,node_len,ksd.depth + node_len
            ));
            continue;
        }

#ifndef DONTMARK
        // make sure this is the only time we see this (k-1)-overlap
        bool already_seen = compare_and_mark_last_k_minus_one_mer (graph.toString(ksd.node), already_extended_from);
        if (already_seen)
        {
            INFO (("... XXX not extending node %s becaues last k-1-mer was already seen\n", seq.toString().c_str()));
            continue;
        }
#endif

        // continue extending from immediately overlapping kmers
        // there may be just one 1 possibility (there was in-branching)

        /** We get the successors of the node. */
        GraphVector<Node> successors = graph.successors (endNode);

        /** We iterate the successors. */
        for (size_t i=0; i<successors.size(); i++)
        {
            kmers_to_traverse.push_back ( NodeDepth<Node> (successors[i].kmer, successors[i].strand, ksd.depth + len_right +1) );
            // ou plutot depth + len_right +1 (+1 = la nt ajoutee ici) (et pas node_len)  ?
        }

        INFO (("... number of extensions: %d\n", successors.size()));

    }   /* end of while (kmers_to_traverse.size() > 0) */

    /** We have to flush the output bank. */
    outputBank->flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Node, typename Edge, typename Graph>
void IterativeExtensions<span, Node, Edge, Graph>::construct_linear_seqs (
    const std::string& L,
    const std::string& R,
    const std::string& output_file,
    bool               swf
)
{
    /** We want to have only one line for data in FASTA output. */
    BankFasta::setDataLineSize (0);

    construct_linear_seqs (L, R, new BankFasta (output_file), swf);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span, typename Node, typename Edge, typename Graph>
void IterativeExtensions<span, Node, Edge, Graph>::buildSequence (
    const Node&     node,
    const Path_t<Node>&     consensusRight,
    size_t          nbNodes,
    size_t          depth,
    Sequence&       seq
)
{
    /** Shortcuts. */
    Data&  data     = seq.getData();
    size_t lenRight = consensusRight.size();

    /** We compute the total size of the sequence. */
    size_t fullLength = graph.getKmerSize() + consensusRight.size();

    /* we need this in mapsembler: the first used kmer should be a k-1 mer. Indeed, if a first kmer is extracted from a sequence :
     * -------------****** (k=6), then this node is the one linked to a new one starting with ******,
     * thus with an overlap of k and not k-1.
     */
    size_t offset = (depth == 0 && dont_output_first_nucleotide) ? 1 : 0;

    /** We compute the total size of the sequence. */
    size_t length = fullLength - offset;

    /** We set the comment of the sequence. */
    seq.setComment (Stringify::format ("%lli__len__%i__depth__%i", nbNodes, fullLength, depth) );

    /** We set the data length. */
    seq.getData().resize (length);

    size_t idx=0;
    size_t i = offset;

    /** We dump the starting node. */
    string nodeStr = graph.toString (node);
    for ( ; i<nodeStr.size(); i++)  { data[idx++] = nodeStr[i]; }

    /** We dump the right part. */
    for (i=0; i<lenRight; i++)  {  data[idx++] = ascii (consensusRight[i]); }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
