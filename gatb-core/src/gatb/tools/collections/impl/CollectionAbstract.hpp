/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Collection.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_ABSTRACT_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_ABSTRACT_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

template <class Item>
class CollectionAbstract : public Collection<Item>
{
public:

    /** Constructor.
     * \param id : identifier of the node
     * \param bag : reference on the bag delegate.
     * \param iterable : reference on the iterable delegate
     */
    CollectionAbstract (Bag<Item>* bag, Iterable<Item>* iterable)
        : _bag(0), _iterable(0)
    {
        setBag      (bag);
        setIterable (iterable);
    }

    /** Destructor. */
    virtual ~CollectionAbstract()
    {
        setBag      (0);
        setIterable (0);
    }

    /** \copydoc Collection::bag */
    Bag<Item>* bag() { return _bag; }

    /** \copydoc Collection::iterable */
    Iterable<Item>* iterable()  { return _iterable; }

    /** \copydoc Iterable::iterator */
    dp::Iterator<Item>* iterator ()  { return _iterable->iterator(); }

    /** \copydoc Iterable::getNbItems */
    int64_t getNbItems ()  { return _iterable->getNbItems(); }

    /** \copydoc Bag::insert */
    void insert (const Item& item)  { _bag->insert (item); }

    /** \copydoc Bag::flush */
    void flush ()  { _bag->flush(); }

protected:

    Bag<Item>* _bag;
    void setBag (Bag<Item>* bag)  { SP_SETATTR(bag); }

    Iterable<Item>* _iterable;
    void setIterable (Iterable<Item>* iterable)  { SP_SETATTR(iterable); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_ABSTRACT_HPP_ */
