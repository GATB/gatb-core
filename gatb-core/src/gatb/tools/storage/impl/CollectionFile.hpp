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

/** \file CollectionFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_FILE_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/CollectionAbstract.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>
#include <vector>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of the Collection interface with a file.
 *
 * This implementation reads/writes Item objects in a file.
 */
template <class Item> class CollectionFile : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionFile (const std::string& filename, size_t cacheItemsNb=10000)
        : collections::impl::CollectionAbstract<Item> (
             new collections::impl::BagFile<Item>(filename),
             new collections::impl::IterableFile<Item>(filename, cacheItemsNb)
          ),  _name(filename)
    {}

    /** Destructor. */
    virtual ~CollectionFile() {}

    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {  gatb::core::system::impl::System::file().remove (_name);  }

private:

    std::string _name;
};

/********************************************************************************/
/* Experimental (not documented). */
template <class Item> class CollectionGzFile : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    CollectionGzFile (const std::string& filename, size_t cacheItemsNb=10000)
    : collections::impl::CollectionAbstract<Item> (
                                                   new collections::impl::BagGzFile<Item>(filename),
                                                   new collections::impl::IterableGzFile<Item>(filename, cacheItemsNb)
                                                   ),  _name(filename)
    {}
    
    /** Destructor. */
    virtual ~CollectionGzFile() {}
    
    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {  gatb::core::system::impl::System::file().remove (_name);  }
    
private:
    
    std::string _name;
};
  
/********************************************************************************/
/* Experimental (not documented). */
template <class Item> class CollectionCountFile : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    CollectionCountFile (const std::string& filename, size_t cacheItemsNb=10000)
    : collections::impl::CollectionAbstract<Item> (
                                                   new collections::impl::BagCountCompressedFile<Item>(filename),
                                                   new collections::impl::IterableCountCompressedFile<Item>(filename, cacheItemsNb)                                                    ),  _name(filename)
    {}
    
    /** Destructor. */
    virtual ~CollectionCountFile() {}
    
    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {  gatb::core::system::impl::System::file().remove (_name);  }
    
private:
    
    std::string _name;
};
    
/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_FILE_HPP_ */
