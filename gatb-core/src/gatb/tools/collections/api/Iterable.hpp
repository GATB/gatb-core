/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
template <class Item> class Iterable //: public dp::SmartPointer
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
