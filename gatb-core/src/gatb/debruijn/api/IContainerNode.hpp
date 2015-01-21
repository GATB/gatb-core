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

/** \file IContainerNode.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Container interface
 *
 *  This file holds interfaces related to the Container interface, ie something we can ask for items.
 */

#ifndef _GATB_CORE_DEBRUIJN_ICONTAINER_NODE_HPP_
#define _GATB_CORE_DEBRUIJN_ICONTAINER_NODE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
/********************************************************************************/

/** \brief Container interface
 *
 * The Container interface provides an operation that ask for a given item.
 *
 * This interface is mainly used by the impl::Graph class
 */
template <class Item> class IContainerNode : public tools::collections::Container<Item>
{
public:

    /** Destructor. */
    virtual ~IContainerNode() {}

    /** Tells whether an item exists or not in the container
     * \return true if the item exists, false otherwise */
    virtual bool contains (const Item& item) = 0;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_ICONTAINER_NODE_HPP_ */
