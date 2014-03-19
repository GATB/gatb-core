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

/** \file ContainerSet.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Container implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_CONTAINER_SET_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_CONTAINER_SET_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/system/api/types.hpp>

#include <vector>
#include <algorithm>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class ContainerSet : public Container<Item>, public system::SmartPointer
{
public:
    
    ContainerSet (dp::Iterator<Item>* it)
    {
        LOCAL (it);
        for (it->first(); !it->isDone(); it->next())  {  _items.push_back (it->item());  }

        std::sort (_items.begin(), _items.end());
    }
    
    /** */
    bool contains (const Item& item)
    {
        return std::binary_search (_items.begin(), _items.end(), item);
    }

private:
    
    std::vector<Item> _items;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_CONTAINER_SET_HPP_ */
