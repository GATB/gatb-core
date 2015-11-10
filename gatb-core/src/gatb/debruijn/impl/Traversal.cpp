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

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <gatb/debruijn/impl/Traversal.hpp>
#include <gatb/debruijn/impl/Frontline.hpp>

using namespace std;
using namespace gatb::core::tools::misc;

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

/********************************************************************************/
#define DEBUG(a)   //a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
TraversalTemplate<Node,Edge,GraphDataVariant>* TraversalTemplate<Node,Edge,GraphDataVariant>::create (
    TraversalKind   type,
    const GraphTemplate<Node,Edge,GraphDataVariant>&    graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&     terminator,
    int             max_len,
    int             max_depth,
    int             max_breadth
)
{
    TraversalTemplate<Node,Edge,GraphDataVariant>* result = 0;

         if (type == TRAVERSAL_UNITIG)  { result = new SimplePathsTraversalTemplate<Node,Edge,GraphDataVariant> (graph, terminator, max_len, max_depth, max_breadth); }
    else if (type == TRAVERSAL_CONTIG)  { result = new MonumentTraversalTemplate<Node,Edge,GraphDataVariant>    (graph, terminator, max_len, max_depth, max_breadth); }
    else if (type == TRAVERSAL_NONE)    { result = new NullTraversalTemplate<Node,Edge,GraphDataVariant>        (graph, terminator, max_len, max_depth, max_breadth); }
    else                                { result = new MonumentTraversalTemplate<Node,Edge,GraphDataVariant>    (graph, terminator, max_len, max_depth, max_breadth); }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
TraversalTemplate<Node,Edge,GraphDataVariant>* TraversalTemplate<Node,Edge,GraphDataVariant>::create (
    const std::string&  type,
    const GraphTemplate<Node,Edge,GraphDataVariant>&        graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>&         terminator,
    int                 max_len,
    int                 max_depth,
    int                 max_breadth
)
{
    TraversalKind typeEnum;  parse (type, typeEnum);

    return create (typeEnum, graph, terminator, max_len, max_depth, max_breadth);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
TraversalTemplate<Node,Edge,GraphDataVariant>::TraversalTemplate (
    const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator,
    int max_len,
    int max_depth,
    int max_breadth
)
    : final_stats(TraversalStats()), stats(TraversalStats()),
      graph(graph), terminator(terminator),
      maxlen      (max_len     == 0 ? TraversalTemplate<Node,Edge,GraphDataVariant>::defaultMaxLen     : max_len),
      max_depth   (max_depth   == 0 ? TraversalTemplate<Node,Edge,GraphDataVariant>::defaultMaxDepth   : max_depth),
      max_breadth (max_breadth == 0 ? TraversalTemplate<Node,Edge,GraphDataVariant>::defaultMaxBreadth : max_breadth)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
int TraversalTemplate<Node,Edge,GraphDataVariant>::traverse (Node& startingNode, Node& currentNode, Direction dir, Path_t<Node>& consensus)
{
    currentNode = startingNode;
    Node previousNode;

    int nnt = 0;

    bool looping = false;

    Path_t<Node> path;  path.resize (this->max_depth+1);

    int bubble_start=0, bubble_end=0;
    bubbles_positions.clear();

    /** We reset the consensus to be filled. */
    consensus.clear ();

    while( (nnt = avance (currentNode, dir, consensus.size() == 0, path, previousNode)))
    {
        // found branching or marked kmers
        if (nnt < 0)  {   break;  }

        if (nnt > 1)  {  bubble_start = consensus.size();  }

        // keep re-walking the nucleotides we just discovered, to append to consensus and mark kmers as we go
        for (int i = 0; i < nnt; i++)
        {
            /** We add the current nucleotide into the contig. */
            consensus.push_back (path[i]);

            /** We update previous and current nodes. */
            previousNode = currentNode;

            /** We compute the neighbor node for the current nucleotide of the path.
             * WARNING: we 'build' here a node without checking that it belongs to the Debruijn graph, because
             * the transition nucleotide has been got through a call to Graph::neighbors, and therefore it should
             * be a trustable transition nucleotide. */
            currentNode  = this->graph.neighbor (currentNode, dir, path[i]);

            /** We mark the node as used in assembly. */
            terminator.mark (currentNode);

            /** perfectly circular regions with no large branching can happen (rarely) */
            if (currentNode.kmer == startingNode.kmer)  {  looping = true;  }
        }

        if (nnt > 1)
        {
            bubble_end = consensus.size();
            bubbles_positions.push_back(std::make_pair(bubble_start,bubble_end));
        }         

        if (looping)  {  break;  }

        if ((int)consensus.size() > maxlen)  {  break;  }

    }  /* end of while( (nnt = avance ... */

    return consensus.size();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
float TraversalTemplate<Node,Edge,GraphDataVariant>::needleman_wunch (const Path_t<Node>& a, const Path_t<Node>& b)
{
    float gap_score = -5;
    float mismatch_score = -5;
    float match_score = 10;
    #define nw_score(x,y) ( (x == y) ? match_score : mismatch_score )

    int n_a = a.size(), n_b = b.size();
    float ** score =  (float **) MALLOC (sizeof(float*) * (n_a+1));
    for (int ii=0; ii<(n_a+1); ii++)
    {
        score [ii] = (float *) MALLOC (sizeof(float) * (n_b+1));
    }


    for (int i = 0; i <= n_a; i++)
        score[i][0] = gap_score * i;
    for (int j = 0; j <= n_b; j++)
        score[0][j] = gap_score * j;

    // compute dp
    for (int i = 1; i <= n_a; i++)
    {
        for (int j = 1; j <= n_b; j++)
        {
            float match = score[i - 1][j - 1] + nw_score(a[i-1],b[j-1]);
            float del =  score[i - 1][j] + gap_score;
            float insert = score[i][j - 1] + gap_score;
            score[i][j] = max( max(match, del), insert);
        }
    }

    // traceback
    int i=n_a, j=n_b;
    float identity = 0;
    while (i > 0 && j > 0)
    {
        float score_current = score[i][j], score_diagonal = score[i-1][j-1], score_up = score[i][j-1], score_left = score[i-1][j];
        if (score_current == score_diagonal + nw_score(a[i-1], b[j-1]))
        {
            if (a[i-1]== b[j-1])
                identity++;
            i -= 1;
            j -= 1;
        }
        else
        {
            if (score_current == score_left + gap_score)
                i -= 1;
            else if (score_current == score_up + gap_score)
                    j -= 1;
        }
    }
    identity /= max( n_a, n_b); // modif GR 27/09/2013    max of two sizes, otherwise free gaps



    for (int ii=0; ii<(n_a+1); ii++)
    {
        FREE (score [ii]);
    }
    FREE(score);

    return identity;
}

/*********************************************************************
            #     #  #     #  ###  #######  ###   #####
            #     #  ##    #   #      #      #   #     #
            #     #  # #   #   #      #      #   #
            #     #  #  #  #   #      #      #   #  ####
            #     #  #   # #   #      #      #   #     #
            #     #  #    ##   #      #      #   #     #
             #####   #     #  ###     #     ###   #####
*********************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
SimplePathsTraversalTemplate<Node,Edge,GraphDataVariant>::SimplePathsTraversalTemplate (
    const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : TraversalTemplate<Node,Edge,GraphDataVariant>(graph, terminator, maxlen, max_depth, max_breadth)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
char SimplePathsTraversalTemplate<Node,Edge,GraphDataVariant>::avance (
    Node& node,
    Direction dir,
    bool first_extension,
    Path_t<Node>& path,
    Node& previousNode
)
{
    int res = this->graph.simplePathAvance (node, dir, path[0]);
    this->deadend = false; // whether the traversal ends with no extension or not -- helps filter out isolated contigs in Minia
    switch (res)
    {
        case -2:
           this->stats.couldnt_inbranching++; break;
        case -1:
           this->stats.couldnt_outbranching++; break;
        case 0:
           this->stats.couldnt_no_extension++; 
           this->deadend = true; break;
    }

    return  max (res,  0);
}

/*********************************************************************
#     #  #######  #     #  #     #  #     #  #######  #     #  #######
##   ##  #     #  ##    #  #     #  ##   ##  #        ##    #     #
# # # #  #     #  # #   #  #     #  # # # #  #        # #   #     #
#  #  #  #     #  #  #  #  #     #  #  #  #  #####    #  #  #     #
#     #  #     #  #   # #  #     #  #     #  #        #   # #     #
#     #  #     #  #    ##  #     #  #     #  #        #    ##     #
#     #  #######  #     #   #####   #     #  #######  #     #     #
*********************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::MonumentTraversalTemplate (
    const GraphTemplate<Node,Edge,GraphDataVariant>& graph,
    TerminatorTemplate<Node,Edge,GraphDataVariant>& terminator,
    int maxlen,
    int max_depth,
    int max_breadth
)
    : TraversalTemplate<Node,Edge,GraphDataVariant>(graph, terminator, maxlen, max_depth, max_breadth)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
char MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::avance (
    Node& node,
    Direction dir,
    bool first_extension,
    Path_t<Node>& consensus,
    Node& previousNode
)
{
    // if we're on a simple path, just traverse it
    int is_simple_path = this->graph.simplePathAvance (node, dir, consensus[0]);
    if (is_simple_path > 0)  {  return 1;  }

    // the following function does:
    // * a bfs from the starting kmer, stopping when:
    //      - breadth > max_breadth
    //      or
    //      - depth > max_depth
    // * check if there a single end point
    // * computing all possible paths between start and end
    // * returns one flattened consensus sequence
    bool success = explore_branching (node, dir, consensus, previousNode);
    if (!success)
    { 
        this->stats.ended_traversals++;
        return 0;
    }

    return consensus.size();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::explore_branching (
    Node& node,
    Direction dir,
    Path_t<Node>& consensus,
    Node& previousNode
)
{
    set<Node> all_involved_extensions;

    return explore_branching (node, dir, consensus, previousNode, all_involved_extensions);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::explore_branching (
    Node& startNode,
    Direction dir,
    Path_t<Node>& consensus,
    Node& previousNode,
    std::set<Node>& all_involved_extensions
)
{
    Node endNode;

    // find end of branching, record all involved extensions (for future marking)
    // it returns false iff it's a complex bubble
    int traversal_depth = find_end_of_branching (dir, startNode, endNode, previousNode, all_involved_extensions);
    if (!traversal_depth)  
    {
        this->stats.couldnt_find_all_consensuses++;
        return false;
    }

    // find all consensuses between start node and end node
    bool success;
    set<Path_t<Node> > consensuses = all_consensuses_between (dir, startNode, endNode, traversal_depth+1, success);

    // if consensus phase failed, stop
    if (!success)  {  return false;  }

    consensus.resize (0);
    // validate paths, based on identity
    bool validated = validate_consensuses (consensuses, consensus);
    if (!validated)   
    {  
        this->stats.couldnt_validate_consensuses++;
        return false;  
    }

    // the consensuses agree, mark all the involved extensions
    // (corresponding to alternative paths we will never traverse again)
    mark_extensions (all_involved_extensions);

    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
//template <typename Frontline=FrontlineBranching> // TODO: someday do something along this line to refactor GraphSimplification in minia
template <typename Node, typename Edge, typename GraphDataVariant>
int MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::find_end_of_branching (
    Direction    dir,
    Node&  startingNode,
    Node&        endNode,
    Node&  previousNode,
    std::set<Node>& all_involved_extensions
)
{
    /** We need a branching frontline. */
    FrontlineBranchingTemplate<Node,Edge,GraphDataVariant> frontline (dir, this->graph, this->terminator, startingNode, previousNode, &all_involved_extensions);

    do  {
        bool should_continue = frontline.go_next_depth();
        if (!should_continue) 
        {
            if (frontline.stopped_reason == FrontlineTemplate<Node,Edge,GraphDataVariant>::MARKED)
                this->stats.couldnt_because_marked_kmer++;
            if (frontline.stopped_reason == FrontlineTemplate<Node,Edge,GraphDataVariant>::IN_BRANCHING_DEPTH)
                this->stats.couldnt_inbranching_depth++;
            if (frontline.stopped_reason == FrontlineTemplate<Node,Edge,GraphDataVariant>::IN_BRANCHING_BREADTH)
                this->stats.couldnt_inbranching_breadth++;
            if (frontline.stopped_reason == FrontlineTemplate<Node,Edge,GraphDataVariant>::IN_BRANCHING_OTHER)
                this->stats.couldnt_inbranching_other++;
            return 0;
        }

        // don't allow a depth too large
        if ((int)frontline.depth() > this->max_depth)
        {  
            this->stats.couldnt_traverse_bubble_depth++;
            return 0;  
        }

        // don't allow a breadth too large
        if ((int)(frontline.size())> this->max_breadth)
        {  
            this->stats.couldnt_traverse_bubble_breadth++;
            return 0;
        }

        // stopping condition: frontline is either empty, or contains only 1 kmer
        // needs the kmer to be non-branching, in order to avoid a special case of bubble immediatly after a bubble
        // affects mismatch rate in ecoli greatly
        if (frontline.size() == 0)  
        {
            this->stats.couldnt_find_extension++;
            return 0;
        }

         if (frontline.size() == 1) // longer contigs but for some reason, higher mismatch rate
             {break;}
            // TODO: didn't implement is_branching in MPHFTerminator, but the line above used to be commented and the one below enabled. see if it still affects assembly.
        //if (frontline.size() == 1 &&   (!terminator.isEnabled() || !terminator.is_branching(frontline.front().node)) )  {   break;  }
    }
    while (1);

    if (frontline.size()==1)
    {
        endNode = frontline.front().node;
        return frontline.depth();
    }

   return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
void MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::mark_extensions (std::set<Node>& extensions_to_mark)
{
    if (this->terminator.isEnabled())
    {
        for(typename std::set<Node>::iterator it = extensions_to_mark.begin(); it != extensions_to_mark.end() ; ++it)
        {
            Node node = *it; // need this, because terminator.mark will want to modify node to cache its mphf index, hence it cannot be const.
            this->terminator.mark (node);
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
set<Path_t<Node> > MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::all_consensuses_between (
    Direction    dir,
    Node& startNode,
    Node& endNode,
    int traversal_depth,
    std::set<typename Node::Value> usedNode,
    Path_t<Node> current_consensus,
    bool& success
)
{
    set<Path_t<Node> > consensuses;

    // find_end_of_branching and all_consensues_between do not always agree on clean bubbles ends
    // until I can fix the problem, here is a fix
    // to reproduce the problem: SRR001665.fasta 21 4
    if (traversal_depth < -1)
    {
        success = false;
        this->stats.couldnt_consensus_negative_depth++;
        return consensuses;
    }

    if (startNode.kmer == endNode.kmer)// not testing for end_strand anymore because find_end_of_branching doesn't care about strands
    {
        consensuses.insert(current_consensus);
        return consensuses;
    }

    /** We retrieve the neighbors of the provided node. */
    typename GraphTemplate<Node,Edge,GraphDataVariant>::template Vector<Edge> neighbors = this->graph.neighborsEdge (startNode, dir);

    /** We loop these neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** Shortcut. */
        Edge& edge = neighbors[i];

        // don't resolve bubbles containing loops
        // (tandem repeats make things more complicated)
        // that's a job for a gapfiller
        if (usedNode.find(edge.to.kmer) != usedNode.end())
        {
            success = false;
            this->stats.couldnt_consensus_loop++;
            return consensuses;
        }

        // generate extended consensus sequence
        Path_t<Node> extended_consensus(current_consensus);
        extended_consensus.push_back (edge.nt);

        // generate list of used kmers (to prevent loops)
        set<typename Node::Value> extended_kmers (usedNode);
        extended_kmers.insert (edge.to.kmer);

        // recursive call to all_consensuses_between
        set<Path_t<Node> > new_consensuses = all_consensuses_between (
            dir,
            edge.to,
            endNode,
            traversal_depth - 1,
            extended_kmers,
            extended_consensus,
            success
        );

        consensuses.insert (new_consensuses.begin(), new_consensuses.end());

        // mark to stop we end up with too many consensuses
        if (consensuses.size() > (unsigned int)this->max_breadth)  {
            this->stats.couldnt_consensus_amount++;
            success = false;  
        }

        // propagate the stop if too many consensuses reached
        if (success == false)  {   return consensuses;  }
    }

    return consensuses;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
set<Path_t<Node> > MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::all_consensuses_between (
    Direction    dir,
    Node& startNode,
    Node& endNode,
    int traversal_depth,
    bool &success
)
{
    set<typename Node::Value> usedNode;
    usedNode.insert(startNode.kmer);
    Path_t<Node> current_consensus;
    current_consensus.start = startNode;
    success = true;

    return all_consensuses_between (dir, startNode, endNode, traversal_depth, usedNode, current_consensus, success);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::validate_consensuses (set<Path_t<Node> >& consensuses, Path_t<Node>& result)
{
    bool debug = false;
    // compute mean and stdev of consensuses
    int mean = 0;
    int path_number = 0;
    for(typename set<Path_t<Node> >::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        //if (debug)  printf("bubble path %d: %s (len=%lu)\n",path_number,(*it).c_str(),(*it).length());
        mean+=(*it).size();
        path_number++;
    }
    mean/=consensuses.size();
    double stdev = 0;
    for(typename set<Path_t<Node> >::iterator it = consensuses.begin(); it != consensuses.end() ; ++it)
    {
        int consensus_length = (*it).size();
        stdev += pow(fabs(consensus_length-mean),2);
    }
    stdev = sqrt(stdev/consensuses.size());

    // don't traverse large bubbles
    if (mean > this->max_depth)
    {
        this->stats.couldnt_validate_bubble_mean_depth++;
        return false;
    }

    // don't traverse large deadends (here, having one consensus means the other paths were large deadends)
    if (consensuses.size() == 1 && mean > (int)(this->graph.getKmerSize())+1) // deadend length should be < k+1 (most have length 1, but have seen up to 10 in ecoli)
    {
        this->stats.couldnt_validate_bubble_deadend++;
        return false;
    }

    if (debug) printf("%lu-bubble mean %d, stdev %.1f\n",consensuses.size(),mean,stdev);

    // traverse bubbles if paths have roughly the same length
    if (stdev>mean/5)
    {
        this->stats.couldnt_validate_bubble_stdev++;
        return false;
    }

    // check that all consensuses are similar
    if (! all_consensuses_almost_identical(consensuses))
    {
        this->stats.couldnt_validate_bubble_identity++;
        return false;
    }

    // if all good, an arbitrary consensus is chosen (if no MPHF) or the most abundance one is chosen (if MPHF available)
    bool has_mphf = this->graph.checkState(GraphTemplate<Node,Edge,GraphDataVariant>::STATE_MPHF_DONE);

    Path_t<Node> chosen_consensus;
    if (has_mphf)
        chosen_consensus = most_abundant_consensus(consensuses);
    else
        chosen_consensus = *consensuses.begin();

    int result_length = chosen_consensus.size();
    if  (result_length> this->max_depth) // it can happen that consensus is longer than max_depth, despite that we didn't explore that far (in a messy bubble with branchings inside)
    {
        this->stats.couldnt_validate_bubble_long_chosen++;
        return false;
    }

    /** We the the result consensus. */
    result = chosen_consensus;
    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
bool MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::all_consensuses_almost_identical (set<Path_t<Node> >& consensuses)
{
    for (typename set<Path_t<Node> >::iterator it_a = consensuses.begin(); it_a != consensuses.end(); it_a++)
    {
        typename set<Path_t<Node> >::iterator it_b = it_a;
        std::advance(it_b,1);
        while (it_b != consensuses.end())
        {
            int identity = this->needleman_wunch(*it_a,*it_b) * 100;
            if (identity < consensuses_identity)
            {
                //cout << "couldn't pop bubble due to identity %:" << identity << " over length " << (*it_a).size() << " " << (*it_b).size() << endl;
                return false;
            }
            std::advance(it_b,1);
        }
    }
    return true;
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : might have a bug, see remark in there. need investigation.
*********************************************************************/
template <typename Node, typename Edge, typename GraphDataVariant>
Path_t<Node> MonumentTraversalTemplate<Node,Edge,GraphDataVariant>::most_abundant_consensus(set<Path_t<Node> >& consensuses)
{
    Path_t<Node> res;
    bool debug = false;

    unsigned long best_mean_abundance = 0;
    string best_p_str = "";

    if (debug)
        cout << endl << "starting to decide which consensus to choose" << endl;

    for (typename set<Path_t<Node> >::iterator it = consensuses.begin(); it != consensuses.end(); it++)
    {
        // iterate over all kmers in consensus and get mean abundance
        Path_t<Node> p = *it;

        // FIXME: I think that code might be buggy!! (wrong p_str constructed in the bubble.fa example of Minia. see GraphSimplification.cpp for a potential fix)
        // it might have to do with a previous DIR_INCOMING nt bug.

        // naive conversion from path to string
        string p_str = this->graph.toString(p.start);
        for (size_t i = 0; i < p.size(); i++)
            p_str.push_back(p.ascii(i));

        if (debug)
            cout << endl << "mean cov for path: " << p_str << endl << "abundance: " << endl;

        unsigned long mean_abundance = 0;
        for (size_t i = 0; i < p.size(); i++)
        {            
            Node node = this->graph.buildNode((char *)(p_str.c_str()), i); 
            /* I know that buildNode was supposed to be used for test purpose only,
             * but couldn't find anything else to transform my substring into a kmer */

            unsigned char abundance = this->graph.queryAbundance(node);
            mean_abundance += abundance;
        
            if (debug)
                cout << (unsigned int)abundance << " ";
        }
        mean_abundance /= p.size();
        
        if (debug)
            cout << "mean: " << mean_abundance << endl;

        if (mean_abundance > best_mean_abundance)
        {
            best_mean_abundance = mean_abundance;
            res = p;
            best_p_str = p_str;
        }

    }

    if (debug && consensuses.size() > 1)
        cout << endl << "chosen path: " << best_p_str << " abundance: " << best_mean_abundance << endl;

    return res;
}

// legacy GATB compatibility
template class TraversalTemplate<Node, Edge, GraphDataVariant>; 
template class MonumentTraversalTemplate<Node, Edge, GraphDataVariant>; 
template class SimplePathsTraversalTemplate<Node, Edge, GraphDataVariant>; 
template class NullTraversalTemplate<Node, Edge, GraphDataVariant>; 


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
