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

/** \file Iterable.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterable interface
 *
 *  This file holds interfaces related to the Iterable interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_ITERABLE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_ITERABLE_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/SmartPointer.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
/** \brief Tools package */
namespace tools         {
/** \brief Collections interfaces */
namespace collections   {
/********************************************************************************/

/** \brief Iterable interface
 *
 * The Iterable interface provides an operation that creates an iterator.
 *
 * Note that one Iterable instance can create several iterators.
 */
template <class Item> class Iterable : public dp::SmartPointer
{
public:

    /** Create an iterator for the given Iterable instance.
     * \return the new iterator. */
    virtual dp::Iterator<Item>* iterator () = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_ITERABLE_HPP_ */
