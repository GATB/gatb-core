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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPHUNITIGS_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPHUNITIGS_HPP_

/********************************************************************************/
#include <vector>
#include <set>

#ifdef USE_NEW_CXX
#include <unordered_map>
#define NS_TR1_PREFIX std
#else
#include <tr1/unordered_map>
#define NS_TR1_PREFIX std::tr1
#endif

#include <gatb/debruijn/impl/Graph.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for De Bruijn graph management. */
namespace debruijn  {
/** \brief Implementation package for De Bruijn graph management. */
namespace impl      {

/********************************************************************************/

    enum Unitig_pos 
    {
        UNITIG_BEGIN = 1,
        UNITIG_END = 2,
        UNITIG_BOTH = 3 /* making sure it's BEGIN | END */
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

/** \brief Class representing a De Bruijn graph on top of unitigs.
 * a bit drastic here: no more big templatization, it's built on top of nodefast and edgefast,
 * and requires just a span template
 */

template <size_t span>
class GraphUnitigsTemplate : public GraphTemplate<NodeFast<span>, EdgeFast<span>, GraphDataVariantFast<span>>
{
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

    /** Build a graph from scratch.
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static GraphUnitigsTemplate  create (tools::misc::IProperties* options)  {  return  GraphUnitigsTemplate (options);  }

    /** Load a graph from some URI.
     * \param[in] uri : the uri to get the graph from
     * \return the loaded graph.
     */
    static GraphUnitigsTemplate  load (const std::string& uri)  {  return  GraphUnitigsTemplate (uri);  }
    
    static tools::misc::IOptionsParser* getOptionsParser (bool includeMandatory=true);

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
    GraphUnitigsTemplate& operator= (const GraphUnitigsTemplate& graph);

    /**********************************************************************/
    /*                     GLOBAL ITERATOR METHODS                        */
    /**********************************************************************/

    /** Creates an iterator over nodes of the graph.
     * \return the nodes iterator. */

    inline GraphIterator<Node> iterator () const  {  return getNodes ();           }
    inline GraphIterator<Node> iteratorCachedNodes () const { return getNodes(); } /* cached nodes are just nodes in this case*/

    /**********************************************************************/
    /*                     ALL NEIGHBORS METHODS                          */
    /**********************************************************************/

    /** Returns a vector of neighbors of the provided node.
     * \param[in] node : the node whose neighbors are wanted
     * \param[in] direction : the direction of the neighbors. If not set, out and in neighbors are computed.
     * \return a vector of the node neighbors (may be empty). 
     * Warning: be sure to check if edge.from (or node.from) is actually your input node, or its reverse complement.
     */
    inline GraphVector<Node> neighbors    ( Node& node, Direction dir=DIR_END) const  {  return getNodes(node, dir);           }
    
    inline Node* neighborsDummy      ( Node& node, Direction dir=DIR_END) const  {   return NULL;           }


    /* neighbors used to be templated.. not anymore, trying to avoid nested template specialization; 
     * so call neighbors for getting nodes, or neighborsEdge for getting edges */
    inline GraphVector<Edge> neighborsEdge    ( Node& node, Direction dir=DIR_END) const  {   return getEdges(node, dir);           }
    inline Edge* neighborsDummyEdge      ( Node& node, Direction dir=DIR_END) const  {   return NULL;           }
    inline GraphVector<Edge> neighborsEdge    ( const typename Node::Value& kmer) const          {  return getEdgeValues (kmer);           }

    /** Shortcut for 'neighbors' method with direction==DIR_OUTCOMING.
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    inline GraphVector<Node> successors   ( Node& node) const                 {  return getNodes(node, DIR_OUTCOMING); }

    /** Shortcut for 'neighbors' method with direction==DIR_INCOMING.
     * \param[in] node : the node whose neighbors are wanted
     * \return a vector of the node neighbors (may be empty).
     */
    inline GraphVector<Node> predecessors ( Node& node) const                 {  return getNodes(node, DIR_INCOMING);  }
    
    inline GraphVector<Edge> successorsEdge   ( Node& node) const                 {  return getEdges(node, DIR_OUTCOMING); }
    inline GraphVector<Edge> predecessorsEdge ( Node& node) const                 {  return getEdges(node, DIR_INCOMING);  }


    /**********************************************************************/
    /*                     ONE NEIGHBOR METHODS                           */
    /**********************************************************************/

    /** Return a specific neighbor from a given node. The neighbor is defined by a direction and the transition
     * nucleotide.
     * IMPORTANT: this method will not check that the neighbor node belongs to the graph: it merely
     * computes the next kmer but doesn't check the Bloom filter. It is supposed that the client has already asked
     * for the neighbors and so knows the valid transitions.
     * \param[in] source : the source neighbor
     * \param[in] dir : the direction of the transition
     * \param[in] nt : the nucleotide of the transition
     * \return the neighbor object.
     */
    inline Node neighbor ( Node& source, Direction dir, kmer::Nucleotide nt) const
    {  bool exists=true; return getNode (source, dir, nt, exists);  }

    /** Return a specific neighbor from a given node. The neighbor is defined by a direction and the transition
     * nucleotide.
     * IMPORTANT: this method will check that the neighbor node belongs to the graph. If the neighbor is not in the
     * graph, the 'exists' parameter is set to false, true otherwise.
     * \param[in] source : the source neighbor
     * \param[in] dir : the direction of the transition
     * \param[in] nt : the nucleotide of the transition
     * \param[out] exists : yes means that the neighbor is in the graph, false otherwise
     * \return the neighbor object.
     */
    inline Node neighbor ( Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
    {  return getNode (source, dir, nt, exists);  }

    /** Shortcut for neighbor with dir==DIR_OUTCOMING. */
    inline Node successor ( Node& source, kmer::Nucleotide nt) const
    {  bool exists=true; return getNode (source, DIR_OUTCOMING, nt, exists);  }

    inline Node successor ( Node& source, kmer::Nucleotide nt, bool& exists) const
    {  return getNode (source, DIR_OUTCOMING, nt, exists);  }

   
    /** Shortcut for neighbor with dir==DIR_INCOMING. */
    inline Node predecessor ( Node& source, kmer::Nucleotide nt) const
    {  bool exists=true; return getNode (source, DIR_INCOMING, nt, exists);  }

    inline Node predecessor ( Node& source, kmer::Nucleotide nt, bool& exists) const
    {  return getNode (source, DIR_INCOMING, nt, exists);  }


    /**********************************************************************/
    /*                      MISC NEIGHBORS METHODS                        */
    /**********************************************************************/

    // would be good to factorize with Graph.hpp but i think it'd require having a GraphAbstract and I'm not ready for that kind of design pattern yet.
    size_t indegree  (Node& node) const;
    size_t outdegree (Node& node) const;
    size_t degree    (Node& node, Direction dir) const;
    void degree    (Node& node, size_t& in, size_t &out) const;
   
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
    int simplePathAvance (Node& node, Direction dir, Edge& output) const;

    /** */
    int simplePathAvance (Node& node, Direction dir) const;

    /** */
    int simplePathAvance (Node& node, Direction dir, kmer::Nucleotide& nt) const;

    /** */
    GraphIterator<Node> simplePath     (Node& node, Direction dir) const  { return getSimpleNodeIterator(node, dir); }
    GraphIterator<Edge> simplePathEdge (Node& node, Direction dir) const  { return getSimpleEdgeIterator(node, dir); }


    // high-level functions that used to be in Simplifications.cpp
    double       simplePathMeanAbundance     (const Node& node, Direction dir) const;
    unsigned int simplePathLength            (const Node& node, Direction dir) const;
    Node         simplePathLastNode          (const Node& node, Direction dir) const;
    void         simplePathDelete            (const Node& node, Direction dir, NodesDeleter<NodeFast<span>, EdgeFast<span>, GraphUnitigsTemplate<span>>& nodesDeleter);

    std::string simplePathSequence (const Node& node, bool& isolatedLeft, bool& isolatedRight) const;


    /**********************************************************************/
    /*                         NODE METHODS                               */
    /**********************************************************************/

    /** Tells whether or not a node belongs to the graph.
     * \param[in] item : the node
     * \return true if the node belongs to the graph, false otherwise. */
    bool contains (const Node& item) const;
                                  
    bool contains (const typename gatb::core::kmer::impl::Kmer<span>::Type& item) const;

    /** Tells whether the provided node is branching or not.
     * \param[in] node : the node to be asked
     * \return true if the node is branching, false otherwise. */
    bool isBranching (Node& node) const;

    /** Return the abundance of a node by querying the perfect hash function 
     * \param[in] node : the node
     * \return the abundance */
    int queryAbundance (Node& node) const;

    /** Return the state of a node by querying the perfect hash function. A node state is either normal, marked, or deleted.
     * \param[in] node : the node or a node index (unsigned long) from the MPHF
     * \return the abundance */
    int queryNodeState (Node& node) const;
    void setNodeState (Node& node, int state) const;
    void resetNodeState () const ;
    void disableNodeState () const ; // see Graph.cpp for explanation

    // deleted nodes, related to NodeState above
    void deleteNode (Node& node) const;
    void deleteNodesByIndex(std::vector<bool> &bitmap, int nbCores = 1, gatb::core::system::ISynchronizer* synchro=NULL) const;
    bool isNodeDeleted(Node& node) const;

    // a direct query to the MPHF
    unsigned long nodeMPHFIndex(Node& node) const;

    void cacheNonSimpleNodes(unsigned int nbCores, bool verbose) ; // dummy for unitigs

    /**********************************************************************/
    /*                         EDGE METHODS                               */
    /**********************************************************************/

    /** Tells whether the provided edge is simple: outdegree(from)==1 and indegree(to)==1
     * \param[in] edge : the edge to be asked
     * \return true if the edge is simple, false otherwise. */
    bool isSimple (Edge& edge) const;

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
        STATE_BCALM2_DONE         = (1<<3),
        STATE_MPHF_DONE           = (1<<6) // to keep compatibility with Traversal and others who check for MPHF, since we support _some_ of the MPHF-like queries (but not on all nodes)
    };
    typedef u_int64_t State; /* this is a global graph state, not to be confused of the state of a node (deleted or not) */
    State _state; // one of the few variables not shared between Graph and GraphUnitigs
    State getState () const { return _state; }
    bool  checkState (StateMask mask) const { return (_state & (State)mask)==(State)mask; }
    State setState   (StateMask mask) { _state |=  mask; return _state; }
    State unsetState (StateMask mask) { _state &= ~mask; return _state; }

    /** Constructor for empty graph.*/
    GraphUnitigsTemplate (size_t kmerSize);

    /** Constructor. Use for GraphUnitigsTemplate creation (ie. DSK + debloom) and filesystem save. */
    GraphUnitigsTemplate (bank::IBank* bank, tools::misc::IProperties* params);

    /** Constructor. Use for GraphUnitigsTemplate creation (ie. DSK + debloom) and filesystem save. */
    GraphUnitigsTemplate (tools::misc::IProperties* params);

    /** Constructor. Use for reading from filesystem. */
    GraphUnitigsTemplate (const std::string& uri);

public: // was private: before, but had many compilation errors during the change from Graph to GraphTemplate. took the easy route, set it to "public:", it solved everything.
   
    /** */
    GraphIterator<Node> getNodes () const;
    
    unsigned char countNeighbors (Node&, Direction) const; // simple and much faster version of getNodes, for degree(), outdegree(), indegree() queries
    void countNeighbors (Node&, size_t&, size_t&) const;  // compute in and out degree at the same time

    /** */
    GraphIterator<Node> getSimpleNodeIterator (Node& node, Direction dir) const;

    /** */
    GraphIterator<Edge> getSimpleEdgeIterator (Node& node, Direction dir) const;

    /** */
    GraphVector<Edge> getEdges (Node source, Direction direction) const;

    /** */
    GraphVector<Node> getNodes (Node &source, Direction direction)  const;

    /** */
    Node getNode (Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const;
    
    GraphVector<Edge> getEdgeValues (const typename Node::Value& kmer) const;


    // core unitigs graph part
    // btw
    // hack so dirty i'd need to call O2 to get it cleaned up
    class ExtremityInfo 
    {
        public:
        uint32_t unitig;
        bool deleted;
        bool rc; // whether the kmer in canonical form appears as rc in the unitig
        Unitig_pos pos; // whether the kmer is at left extremity of unitig or right extremity
        ExtremityInfo() {} // needed so that the type can be used as values in a hash table
        ExtremityInfo(uint32_t u, bool d, bool r, Unitig_pos p) : unitig(u),deleted(d),rc(r), pos(p) {}
        std::string toString() const
        { return " rc:" + std::to_string(rc) + " p:" + ((pos&UNITIG_BEGIN)?"left":"") + ((pos&UNITIG_END)?"right":"") + " " + " d:" + std::to_string(deleted); }
    };
    
    typedef typename gatb::core::kmer::impl::Kmer<span>::Type           Type;

    // structure that links each kmer to an unitig
    // also used to enumerate kmers
    typedef typename NS_TR1_PREFIX::unordered_map<Type, ExtremityInfo> NodeMap;


    void load_unitigs(std::string unitigs_filename);

    bool node_in_same_orientation_as_in_unitig(const Node& node, const ExtremityInfo& e) const;
    
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::ModelDirect ModelDirect;

    // don't forget to copy those variable sin operator= !!
    Model       *modelK;
    ModelDirect *modelKdirect;
    NodeMap utigs_map;
    std::vector<std::string> unitigs;
};
  

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
