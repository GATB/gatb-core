/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
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

#include <gatb/tools/designpattern/api/SmartPointer.hpp>
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
template <class Item> class Bag : public dp::SmartPointer
{
public:

    /** Insert an item into the bag.
     * \param[in] item : the item to be inserted. */
    virtual void insert (const Item& item) = 0;

    /** Flush the current content. May be useful for implementation that uses a cache. */
    virtual void flush () = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_ */
