/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
