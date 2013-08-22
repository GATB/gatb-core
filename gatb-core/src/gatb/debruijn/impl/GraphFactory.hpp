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

#include <gatb/debruijn/api/IGraph.hpp>
#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/GraphBasic.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/misc/impl/Property.hpp>

#include <stdarg.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {

/********************************************************************************/

class GraphFactory
{
public:

    template<typename T>
    static Graph<T>  createGraph (tools::misc::IProperty* prop, ...)
    {
        tools::misc::IProperties* props = new tools::misc::impl::Properties ();
        LOCAL (props);

        va_list args;
        va_start (args, prop);
        props->add (prop, args);
        va_end (args);

        /** */
        return  Graph<T> (singleton().newGraph<T> (props));
    }

    template<typename T>
    static Graph<T>  createGraph (
        tools::collections::Iterable<T>* solidKmers,
        tools::collections::Iterable<T>* cFPKmers,
        size_t kmerSize
    )
    {
        return  Graph<T> (new GraphBasic<T> (solidKmers, cFPKmers, kmerSize));
    }

    template<typename T>
    static Graph<T>  createGraph (
        bank::IBank* bank,
        size_t       kmerSize,
        size_t       nks=1
    )
    {
        return  Graph<T> (new GraphBasic<T> (bank, kmerSize, nks));
    }


protected:

    static GraphFactory& singleton()  { static GraphFactory instance;  return instance; }

    template<typename T>
    Graph<T>  newGraph (tools::misc::IProperties* props)
    {
        /** */
        return  new GraphBasic<T> (props);
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_GRAPH_FACTORY_HPP_ */
