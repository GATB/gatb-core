/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file GraphBasic.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPH_HPP_

/********************************************************************************/

#include <gatb/debruijn/api/IGraph.hpp>

#include <gatb/bank/api/IBank.hpp>

#include <gatb/tools/collections/impl/ProductFile.hpp>
#include <gatb/tools/collections/impl/ProductHDF5.hpp>

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {

/********************************************************************************/

#if 0
#define ProductFactoryLocal tools::collections::impl::ProductFileFactory
#else
#define ProductFactoryLocal tools::collections::impl::ProductHDF5Factory
#endif

/********************************************************************************/

class Graph : public IGraph
{
public:

    /** Destructor. */
    ~Graph ();

    /** */
    void remove ();

    /** From Container interface. */
    bool contains (const Node& item) const;

    /** Creates an iterator over all nodes of the graph.
     * \return the all nodes iterator. */
    INodeIterator* nodes () const;

    /** */
    size_t getOutEdges (const Node& node, EdgeSet& edges)  const   { return getEdges (node, edges, DIR_OUTCOMING); }

    /** */
    size_t getInEdges (const Node& node, EdgeSet& edges)   const   { return getEdges (node, edges, DIR_INCOMING); }

    /** */
    size_t getSuccessors (const Node& node, NodeSet& nodes) const  { return getNodes (node, nodes, DIR_OUTCOMING); }

    /** */
    size_t getPredecessors (const Node& node, NodeSet& nodes) const { return getNodes (node, nodes, DIR_INCOMING); }

    /** */
    bool isEdge (const Node& u, const Node& v) const { return false; }

    /** */
    tools::misc::IProperties& getInfo () const { return (tools::misc::IProperties&)_info; }

    /** */
    std::string toString (const Node& node, kmer::Strand strand = kmer::STRAND_ALL, int mode=0) const;

    /** */
    void getNearestBranchingRange (const Node& node, Node& begin, Node& end) const;

private:

    /** Constructor. Use for Graph creation (ie. DSK + debloom) and filesystem save. */
    Graph (bank::IBank* bank, tools::misc::IProperties* params);

    /** Constructor. Use for Graph creation (ie. DSK + debloom) and filesystem save. */
    Graph (tools::misc::IProperties* params);

    /** Constructor. Use for reading from filesystem. */
    Graph (const std::string& uri);

    /** Product. */
    tools::collections::impl::Product<ProductFactoryLocal>* _product;
    void setProduct (tools::collections::impl::Product<ProductFactoryLocal>* product)  { SP_SETATTR(product); }
    tools::collections::impl::Group<ProductFactoryLocal>& getProduct(const std::string name="")  { return (*_product) (name); }

    /** Creation information. */
    tools::misc::impl::Properties _info;

    /** Defined as a void* for hiding implementation in cpp file. */
    void* _variant;

    /** */
    size_t getEdges (const Node& source, EdgeSet& edges, Direction direction) const;

    /** */
    size_t getNodes (const Node& source, NodeSet& nodes, Direction direction) const;

    /** */
    void executeAlgorithm (tools::misc::impl::Algorithm& algorithm, tools::misc::IProperties* props, tools::misc::IProperties& info);

    /** Friends. */
    template<typename T> friend class GraphFactoryImpl;
    friend class GraphFactory;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_BASIC_HPP_ */
