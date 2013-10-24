/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file GraphFactory.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPH_FACTORY_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPH_FACTORY_HPP_

/********************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/bank/api/IBank.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {

/********************************************************************************/
/** \brief Factory that creates Graph objects
 *
 * Users have to use this factory for getting Graph instances.
 *
 * Two ways:
 *      1) create a graph from scratch => launch DSK + Debloom and save the result in filesystem
 *      2) load a graph from filesystem
 */
class GraphFactory
{
public:

    /** Build a graph from a given bank.
     * \param[in] bank : bank to get the reads from
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static Graph  create (bank::IBank* bank, tools::misc::IProperties* options)  {  return  Graph (bank, options);  }

    /** Build a graph from scratch.
     * \param[in] options : user parameters for building the graph.
     * \return the created graph.
     */
    static Graph  create (tools::misc::IProperties* options)  {  return  Graph (options);  }

    /** Load a graph from some URI.
     * \parm[in] uri : the uri to get the graph from
     * \return the loaded graph.
     */
    static Graph  load (const std::string& uri)  {  return  Graph (uri);  }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_FACTORY_HPP_ */
