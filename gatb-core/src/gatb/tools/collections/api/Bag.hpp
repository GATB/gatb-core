/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Bag.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bag interface
 *
 *  This file holds interfaces related to the Bag interface, ie something we can put items into.
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/********************************************************************************/

/** \brief Bag interface
 *
 * The Iterable interface provides an operation that creates an iterator.
 *
 * Note that one Iterable instance can create several iterators.
 */
template <class Item> class Bag : public virtual  system::ISmartPointer
{
public:

    /** Insert an item into the bag.
     * \param[in] item : the item to be inserted. */
    virtual void insert (const Item& item) = 0;


    /** Insert items into the bag.
     * \param[in] items : items to be inserted.
     * \param[in] length : the number of items to be inserted. If 0 (default value),
     *                     all the items of the vector are inserted. */
    virtual void insert (const std::vector<Item>& items, size_t length=0)
    {
        size_t n = length==0 ? items.size() : length;
        for (size_t i=0; i<n; i++)  {  insert (items[i]); }
    }

    /** Insert items into the bag.
     * \param[in] items : items to be inserted.
     * \param[in] length : number of items to be inserted. */
    virtual void insert (const Item* items, size_t length)
    {
        for (size_t i=0; i<length; i++)  {  insert (items[i]); }
    }

    /** Flush the current content. May be useful for implementation that uses a cache. */
    virtual void flush () = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_ */
