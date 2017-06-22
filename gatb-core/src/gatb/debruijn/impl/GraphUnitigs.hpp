/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014-2016 
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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPHUNITIGS_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPHUNITIGS_HPP_

/********************************************************************************/
#include <vector>
#include <set>

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/UnitigsConstructionAlgorithm.hpp>
#include <gatb/debruijn/impl/ExtremityInfo.hpp>

#include <gatb/debruijn/impl/dag_vector.hpp> // TODO move it to 3rd party


/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for De Bruijn graph management. */
namespace debruijn  {
/** \brief Implementation package for De Bruijn graph management. */
namespace impl      {

/********************************************************************************/

/* Nodes. Those correspond to left or right extremities of a unitig
 *
 * Big difference with Graph.hpp: in GraphUnitig, we don't store kmers in nodes. Kmers are inferred from unitigs
 */
struct NodeGU
{
    /** Default constructor. */
    NodeGU() : unitig(0), pos(UNITIG_BEGIN), strand(kmer::STRAND_FORWARD)  {}

    /** Constructor.
     */
    NodeGU (const uint64_t unitig, Unitig_pos pos, kmer::Strand strand)
        : unitig(unitig), pos(pos), strand(strand) {}

    NodeGU (const uint64_t unitig, Unitig_pos pos)
        : unitig(unitig), pos(pos), strand(kmer::STRAND_FORWARD) {}


    uint64_t unitig;
    Unitig_pos pos;

    /** Strand telling how to interpret the node in the bi-directed DB graph. */
    kmer::Strand strand;

    /** Overload of operator ==  NOTE: it doesn't care about the strand!!! */
    bool operator== (const NodeGU& other) const  { return unitig == other.unitig && pos == other.pos; }

    bool operator!= (const NodeGU& other) const  { return unitig != other.unitig || pos != other.pos; }

    // this need to be implemented, for traversedNodes.find() in Simplifications
    bool operator< (const NodeGU& other) const  { return (unitig < other.unitig || (unitig == other.unitig && pos < other.pos)); }

    void set (uint64_t unitig, Unitig_pos pos, kmer::Strand strand)
    {
        this->unitig = unitig;
        this->strand = strand;
        this->pos    = pos;
    }

    void reverse()
    {
        strand = StrandReverse(strand);
    }

};

struct EdgeGU 
{
    /** The source node of the edge. */
    NodeGU from;
    /** The target node of the edge. */
    NodeGU to;
    /** The direction of the transition. */
    Direction        direction;
    
    // this need to be implemented, for something in in Simplifications
    bool operator< (const EdgeGU& other) const  { return ((from < other.from) || (from == other.from && to < other.to)); } 

    /** Setter for some attributes of the Edge object.
     * \param[in] unitig_from
     * \param[in] pos_from
     * \param[in] strand_from : strand of the 'from' Node
     * \param[in] unitig_to 
     * \param[in] pos_to
     * \param[in] strand_to : strand of the 'from' Node
     * \param[in] dir : direction of the transition.
     */
    void set (
        uint64_t unitig_from, Unitig_pos pos_from, kmer::Strand strand_from,
        uint64_t unitig_to,   Unitig_pos pos_to,   kmer::Strand strand_to,
        Direction dir
    )
    {
        from.set (unitig_from, pos_from, strand_from);
        to.set   (unitig_to,   pos_to,   strand_to);
        direction = dir;
    }
};

/********************************************************************************
                 #####   ######      #     ######   #     #
                #     #  #     #    # #    #     #  #     #
                #        #     #   #   #   #     #  #     #
                #  ####  ######   #     #  ######   #######
                #     #  #   #    #######  #        #     #
                #     #  #    #   #     #  #        #     #
                 #####   #     #  #     #  #        #     #
********************************************************************************/

/** \brief Class representing a De Bruijn graph based on unitigs.
 *
 * Nodes are k-mer extremities of unitigs. They're represented implicitly internally.
 * Edges are (k-1) overlaps between nodes in different unitigs
 */

template <size_t span>
class GraphUnitigsTemplate : public GraphTemplate<NodeFast<span>, EdgeFast<span>, GraphDataVariantFast<span>>
{
    // but actually.. we almost don't use those nodes now! we use NodeGU.
        typedef NodeFast<span> Node;
        typedef EdgeFast<span> Edge;
        
public:

    /********************************************************************************/
    /*                            STATIC METHODS   (create/load)                    */
    /********************************************************************************/

    /** Build an empty graph.
     * \param[in] kmerSize: kmer size
     * \return the created graph.
     */
    static GraphUnitigsTemplate  create (size_t kmerSize)  {  return  GraphUnitigsTemplate(kmerSize);  }

    /** Build a graph from a given bank.
     * \param[in] bank : bank to get the reads from
     * \param[in] fmt : printf-like format for the command line string
     * \return the created graph.
     */
    static GraphUnitigsTemplate  create (bank::IBank* bank, const char* fmt, ...);

    /** Build a graph from user options.
     * \param[in] fmt: printf-like format
     * \return the created graph.
     */
    static GraphUnitigsTemplate  create (const char* fmt, ...);

    /** so, hm, what's the point of a create() function that just calls a constructor? I really don't get the factory pattern yet
     */
    static GraphUnitigsTemplate  create (tools::misc::IProperties* options, bool load_unitigs_after = true /* will be set to false by BCALM 2*/)  {  return  GraphUnitigsTemplate (options, load_unitigs_after);  }

    /** Load a graph from some URI.
     * \param[in] uri : the uri to get the graph from
     * \return the loaded graph.
     */
    static GraphUnitigsTemplate  load (const std::string& uri)  {  return  GraphUnitigsTemplate (uri);  }
    
    // we don't provide an option parser. use Graph's one
    //static tools::misc::IOptionsParser* getOptionsParser (bool includeMandatory=true);

    /********************************************************************************/
    /*                               CONSTRUCTORS                                   */
    /********************************************************************************/

    /* Default Constructor.*/
    GraphUnitigsTemplate ();

    /* Copy Constructor.*/
    GraphUnitigsTemplate (const GraphUnitigsTemplate& graph);

    /* Destructor. */
    ~GraphUnitigsTemplate ();

    /** Affectation overload. */
    GraphUnitigsTemplate& operator= (GraphUnitigsTemplate const& graph);
    GraphUnitigsTemplate& operator= (GraphUnitigsTemplate&& graph);

    /**********************************************************************/
    /*                     GLOBAL ITERATOR METHODS                        */
    /**********************************************************************/

    /** Creates an iterator over nodes of the graph.
     * \return the nodes iterator. */

    inline GraphIterator<NodeGU> iterator () const  {  return getNodes ();           }
    inline GraphIterator<NodeGU> iteratorCachedNodes () const { return getNodes(); } /* cached nodes are just nodes in this case*/

    /**********************************************************************/
    /*                     ALL NEIGHBORS METHODS                          */
    /**********************************************************************/

    /** Returns a vector of neighbors of the provided node.
     * \param[in] node : the node whose neighbors are wanted
     * \param[in] direction : the direction of the neighbors. If not set, out and in neighbors are computed.
     * \return a vector of the node neighbors (may be empty). 
     * Warning: be sure to check if edge.from (or node.from) is actually your input node, or its reverse complement.
     */
    inline GraphVector<NodeGU> neighbors    ( NodeGU& node, Direction dir=DIR_END) const  {  return getNodes(node, dir);           }
    
    inline NodeGU* neighborsDummy      ( NodeGU& node, Direction dir=DIR_END) const  {   return NULL;           }


    /* neighbors used to be templated.. not anymore, trying to avoid nested template specialization; 
     * so call neighbors for getting nodes, or neighborsEdge for getting edges */
    inline GraphVector<EdgeGU> neighborsEdge    ( NodeGU& node, Direction dir=DIR_END) const  {   return getEdges(node, dir);           }
    inline EdgeGU* neighborsDummyEdge      ( NodeGU& node, Direction dir=DIR_END) const  {   return NULL;           }

    /** Shortcut for 'neighbors' method
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    inline GraphVector<NodeGU> successors   ( NodeGU& node) const                 {  return getNodes(node, DIR_OUTCOMING); }
    inline GraphVector<NodeGU> predecessors ( NodeGU& node) const                 {  return getNodes(node, DIR_INCOMING);  }
    inline GraphVector<EdgeGU> successorsEdge   ( NodeGU& node) const                 {  return getEdges(node, DIR_OUTCOMING); }
    inline GraphVector<EdgeGU> predecessorsEdge ( NodeGU& node) const                 {  return getEdges(node, DIR_INCOMING);  }

    /**********************************************************************/
    /*                      MISC NEIGHBORS METHODS                        */
    /**********************************************************************/

    // would be good to factorize with Graph.hpp but i think it'd require having a GraphAbstract and I'm not ready for that kind of design pattern yet.
    size_t indegree  (const NodeGU& node) const;
    size_t outdegree (const NodeGU& node) const;
    size_t degree    (const NodeGU& node, Direction dir) const;
    void degree      (const NodeGU& node, size_t& in, size_t &out) const;
   
    /**********************************************************************/
    /*                      SIMPLIFICATION METHODS                        */
    /**********************************************************************/

    /* perform tip removal, bulge removal and EC removal, as in Minia */
    void simplify(unsigned int nbCores = 1, bool verbose=true);

    /**********************************************************************/
    /*                         SIMPLE PATH METHODS                        */
    /**********************************************************************/

    /** Simple paths traversal
     *  invariant: the input kmer has no in-branching.
     * \returns
     *       1 if a good extension is found
     *       0 if a deadend was reached
     *      -1 if out-branching was detected
     *      -2 if no out-branching but next kmer has in-branching
     */
    int simplePathAvance (const NodeGU& node, Direction dir, EdgeGU& output) const;
    int simplePathAvance (const NodeGU& node, Direction dir) const;

    /** */
    GraphIterator<NodeGU> simplePath     (const NodeGU& node, Direction dir) const  { return getSimpleNodeIterator(node, dir); }
    GraphIterator<EdgeGU> simplePathEdge (const NodeGU& node, Direction dir) const  { return getSimpleEdgeIterator(node, dir); }


    /**********************************************************************/
    /*                         UNITIGS METHODS (the essential stuff)      */
    /**********************************************************************/

    // convention: unitigXXX works on the unitigs as computed as bcalm. never leaves that unitig
    //             simplepathXXX may traverse multiple unitigs
    bool isLastNode                          (const NodeGU& node, Direction dir) const;
    bool isFirstNode                         (const NodeGU& node, Direction dir) const;
    NodeGU             unitigLastNode          (const NodeGU& node, Direction dir) const;
    NodeGU         simplePathLastNode          (const NodeGU& node, Direction dir) ; /* cannot be const becuse it called Longuest_avance that is sometimes not const.. grr. */
    unsigned int     unitigLength            (const NodeGU& node, Direction dir) const;
    unsigned int simplePathLength            (const NodeGU& node, Direction dir) ; /* same reason as above*/ /* NOTE: return number of traversed kmers, so (nucleotide length-k)*/
    double           unitigMeanAbundance     (const NodeGU& node) const;
    double       simplePathMeanAbundance     (const NodeGU& node, Direction dir) ;
    void             unitigDelete          (NodeGU& node, Direction dir, NodesDeleter<NodeGU, EdgeGU, GraphUnitigsTemplate<span>>& nodesDeleter);
    void             unitigDelete          (NodeGU& node) ;
    void         simplePathDelete          (NodeGU& node, Direction dir, NodesDeleter<NodeGU, EdgeGU, GraphUnitigsTemplate<span>>& nodesDeleter);
    std::string  unitigSequence            (const NodeGU& node, bool& isolatedLeft, bool& isolatedRight) const;
    void         unitigMark                (const NodeGU& node); // used to flag simple path as traversed, in minia
    bool         unitigIsMarked        (const NodeGU& node) const;
    
    std::string simplePathBothDirections(const NodeGU& node, bool& isolatedLeft, bool& isolatedRight, bool dummy, float& coverage);
    // aux function, not meant to be called from outside, but maybe it could.
    void simplePathLongest_avance(const NodeGU& node, Direction dir, int& seqLength, int& endDegree, bool markDuringTraversal, float& coverage, std::string* seq = nullptr, std::vector<NodeGU> *unitigNodes = nullptr) ; 

    void debugPrintAllUnitigs() const;

    NodeGU debugBuildNode(std::string startKmer) const;

    /**********************************************************************/
    /*                         NODE METHODS                               */
    /**********************************************************************/

    std::string toString(const NodeGU& node) const;

    // those are not implemented but I need them headers for compatibility with original Graph
    bool contains (const NodeGU& item) const;
    bool isBranching (const NodeGU& node) const;
    int queryAbundance (const NodeGU& node) const;
    int queryNodeState (const NodeGU& node) const;
    void setNodeState (const NodeGU& node, int state) const;
    void resetNodeState () const ;
    void disableNodeState () const ;
    void deleteNodesByIndex(std::vector<bool> &bitmap, int nbCores = 1, gatb::core::system::ISynchronizer* synchro=NULL) const;
    unsigned long nodeMPHFIndex(const NodeGU& node) const;
    void cacheNonSimpleNodes(unsigned int nbCores, bool verbose); 


    // deleted nodes, related to NodeState above
    void deleteNode (/* cannot be const because nodeDeleter isn't */ NodeGU& node) ;
    bool isNodeDeleted(const NodeGU& node) const;

    /**********************************************************************/
    /*                         EDGE METHODS                               */
    /**********************************************************************/

    /** Tells whether the provided edge is simple: outdegree(from)==1 and indegree(to)==1
     * \param[in] edge : the edge to be asked
     * \return true if the edge is simple, false otherwise. */
    bool isSimple (EdgeGU& edge) const;

    /**********************************************************************/
    /*                         MISC METHODS                               */
    /**********************************************************************/

    /** Remove physically a graph. */
    void remove ();

    /**********************************************************************/
    /*                         TYPES                                      */
    /**********************************************************************/
    enum StateMask
    {
        STATE_INIT_DONE           = (1<<0),
        STATE_CONFIGURATION_DONE  = (1<<1),
        STATE_SORTING_COUNT_DONE  = (1<<2),
        STATE_MPHF_DONE           = (1<<6), // to keep compatibility with Traversal and others who check for MPHF, since we support _some_ of the MPHF-like queries (but not on all nodes)
        STATE_BCALM2_DONE         = (1<<20)
    };
    typedef u_int64_t State; /* this is a global graph state, not to be confused of the state of a node (deleted or not) */
    typedef GraphTemplate<NodeFast<span>, EdgeFast<span>, GraphDataVariantFast<span>> BaseGraphForState;
    State getState () const { return BaseGraphForState::_state; }
    bool  checkState (StateMask mask) const { return (BaseGraphForState::_state & (State)mask)==(State)mask; }
    State setState   (StateMask mask) { BaseGraphForState::_state |=  mask; return BaseGraphForState::_state; }
    State unsetState (StateMask mask) { BaseGraphForState::_state &= ~mask; return BaseGraphForState::_state; }

    /** Constructor for empty graph.*/
    GraphUnitigsTemplate (size_t kmerSize);

    GraphUnitigsTemplate (bank::IBank* bank, tools::misc::IProperties* params);

    GraphUnitigsTemplate (tools::misc::IProperties* params, bool load_unitigs_after);

    /** Constructor. Use for reading from filesystem. */
    GraphUnitigsTemplate (const std::string& uri);

public: // was private: before, but had many compilation errors during the change from Graph to GraphTemplate. took the easy route, set it to "public:", it solved everything.
   
    /** */
    GraphIterator<NodeGU> getNodes () const;
    
    unsigned char countNeighbors (const NodeGU&, Direction) const; // simple and much faster version of getNodes, for degree(), outdegree(), indegree() queries
    void countNeighbors (const NodeGU&, size_t&, size_t&) const;  // compute in and out degree at the same time

    /** */
    GraphIterator<NodeGU> getSimpleNodeIterator (const NodeGU& node, Direction dir) const;

    /** */
    GraphIterator<EdgeGU> getSimpleEdgeIterator (const NodeGU& node, Direction dir) const;

    /** */
    GraphVector<EdgeGU> getEdges (const NodeGU& source, Direction direction) const;

    /** */
    GraphVector<NodeGU> getNodes (const NodeGU &source, Direction direction)  const;

    /** */
    NodeGU getNode (const NodeGU& source, Direction dir, kmer::Nucleotide nt, bool& exists) const;
    
    typedef typename gatb::core::kmer::impl::Kmer<span>::Type           Type;

    void build_unitigs_postsolid(std::string unitigs_filename, tools::misc::IProperties* props);
    void load_unitigs(std::string unitigs_filename);

    void load_unitigs_from_gfa(std::string gfa_filename, unsigned int& kmerSize);
    void print_unitigs_mem_stats(uint64_t avg_incoming_size, uint64_t avg_outcoming_size, uint64_t total_unitigs_size, uint64_t nb_utigs_nucl = 0, uint64_t nb_utigs_nucl_mem = 0);

    bool node_in_same_orientation_as_in_unitig(const NodeGU& node) const;
      
    // support for 2-bit compression of unitigs
    std::string internal_get_unitig_sequence(unsigned int unitig_id) const;
    unsigned int internal_get_unitig_length(unsigned int unitig_id) const;
    std::string internal_compress_unitig(std::string seq) const;

    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::ModelDirect ModelDirect;
 
    // !!!!
    // all member variables should be below this point, because i want you to also read this:
    // don't forget to copy those variables in operator= (and the move operator) 
    // classic source of bugs but i couldn't find a foolproof way.
    // !!!!
    std::vector<uint64_t> incoming, outcoming, incoming_map, outcoming_map;
    dag::dag_vector dag_incoming, dag_outcoming, dag_incoming_map, dag_outcoming_map;
    std::vector<std::string> unitigs;
    std::string packed_unitigs;
    std::vector<uint32_t> unitigs_sizes; 
    dag::dag_vector packed_unitigs_sizes;
    std::vector<float> unitigs_mean_abundance;
    //dag::dag_vector unitigs_sizes;// perf hit: from 45s to 74s in chr14; that's because unitigs_sizes is queried _a lot_ just to check if a unitig is just of length k. could save that space with a bit vector, and actually, just use packed_unitigs_sizes for the rest. so.. just to keep in mind that this is a "todo opt" in case we really want to save the space of unitigs_sizes
    //dag::dag_vector unitigs_mean_abundance; // not a big gain and different assembly quality, so i'm keeping it as vector<float>
    std::vector<bool> unitigs_deleted; // could also be replaced by modifying incoming and outcoming vectors. careful not to affect the prefix sum scheme tho.
    std::vector<bool> unitigs_traversed;
    uint64_t nb_unitigs, nb_unitigs_extremities;
    bool compress_navigational_vectors;
    bool pack_unitigs;
    // !!!!
    // read above
    // !!!!
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
