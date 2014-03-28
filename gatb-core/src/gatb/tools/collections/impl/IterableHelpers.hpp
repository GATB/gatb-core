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

/** \file IterableHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterator implementation for file
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

template<class Type, class Listener=tools::misc::impl::ProgressTimer>
class ProgressIterator : public tools::dp::impl::SubjectIterator<Type>
{
public:
    ProgressIterator (Iterable<Type>& iterable, const char* msg = "progress", size_t divide=100)
        : tools::dp::impl::SubjectIterator<Type> (
            iterable.iterator(),
            (iterable.getNbItems() >= 0 ? iterable.getNbItems() : iterable.estimateNbItems()) / divide,
            new Listener ((iterable.getNbItems() >= 0 ? iterable.getNbItems() : iterable.estimateNbItems()), msg)
    ) {}
};


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_ */
