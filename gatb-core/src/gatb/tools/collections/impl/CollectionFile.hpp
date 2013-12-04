/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file CollectionFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_FILE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/CollectionAbstract.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>
#include <vector>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** \brief Collection interface
 */
template <class Item> class CollectionFile : public CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionFile (const std::string& filename, size_t cacheItemsNb=10000)
        : CollectionAbstract<Item> (
             new BagFile<Item>(filename),
             new IterableFile<Item>(filename, cacheItemsNb)
          ),  _name(filename)
    {}

    /** Destructor. */
    virtual ~CollectionFile() {}

    /** \copydoc Collection::remove */
    void remove ()  {  gatb::core::system::impl::System::file().remove (_name);  }

private:

    std::string _name;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_FILE_HPP_ */
