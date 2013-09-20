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

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_COLLECTION_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_COLLECTION_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/api/Bag.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/********************************************************************************/

/** \brief Collection interface
 *
 * The Collection interface is the union of a Bag and an Iterable interfaces
 */
template <class Item> class Collection : public Bag<Item>, public Iterable<Item>
{
public:

    /** Destructor. */
    virtual ~Collection () {}

    /** \return the bag instance. */
    virtual Bag<Item>* bag() = 0;

    /** \return the iterable instance. */
    virtual Iterable<Item>* iterable() = 0;

    /** */
    virtual void remove () = 0;

    /** */
    virtual void addProperty (const std::string& key, const std::string value) = 0;

    /** */
    virtual std::string getProperty (const std::string& key) = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_COLLECTION_HPP_ */
