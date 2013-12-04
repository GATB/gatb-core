/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file CollectionCache.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_CACHE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_CACHE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/system/impl/System.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** \brief Collection interface
 *
 * The Collection interface is the union of a Bag and an Iterable interfaces
 */
template <class Item> class CollectionCache : public CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionCache (Collection<Item>& ref,  size_t cacheSize, system::ISynchronizer* synchro)
        : CollectionAbstract<Item> (
            new BagCache<Item> (ref.bag(), cacheSize, synchro),
            ref.iterable()
        ), _ref(ref)  {}

    /** Destructor. */
    virtual ~CollectionCache() {}

    /** */
    void remove ()  { _ref.remove(); }

    /** */
    Collection<Item>& getRef ()  { return _ref; }

private:
    Collection<Item>& _ref;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_CACHE_HPP_ */
