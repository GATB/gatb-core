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

    ProgressIterator (tools::dp::Iterator<Type>* iterator, const char* msg, size_t nbItems)
        : tools::dp::impl::SubjectIterator<Type> (
            iterator,
            nbItems / 100,
            new Listener (nbItems, msg)
    ) {}
};

/********************************************************************************/

/** \brief Adaptor of an Iterable<T1> into an Iterable<T2>
 */
template <class T1, class T2, class Adaptor>
class IterableAdaptor : public Iterable<T2>, public system::SmartPointer
{
public:
    /** */
    IterableAdaptor (Iterable<T1>& ref)  : _ref(ref) {}

    /** Create an iterator for the given Iterable instance.
     * \return the new iterator. */
    dp::Iterator<T2>* iterator ()  { return new tools::dp::impl::IteratorAdaptor<T1,T2,Adaptor> (_ref.iterator()); }

    /** Return the number of items. If a specific implementation doesn't know the value,
     * it should return -1 by convention.
     * \return the number of items if known, -1 otherwise. */
    int64_t getNbItems () { return _ref.getNbItems(); }

    /** Return the (approximate) number of items. If a specific implementation doesn't know the value,
     * it should return -1 by convention.
     * \return the number of items if known, -1 otherwise. */
    int64_t estimateNbItems ()  { return _ref.estimateNbItems(); }

    /** Return a buffer of items.
     * \param[out] buffer : the buffer
     * \return the buffer */
    T2* getItems (T2*& buffer) { return 0; } //_ref.getItems (buffer); }

    /** */
    size_t getItems (T2*& buffer, size_t start, size_t nb) { return 0; }//_ref.getItems (buffer, start, nb); }

private:
    Iterable<T1>& _ref;
};

/********************************************************************************/

/** \brief some helper methods for Iterable objects. */
class IterableHelpers
{
public:

    template<typename T>
    static bool getItems (Iterable<T>& iterable, std::vector<T>& items)
    {
        size_t nbItems = items.size();

        dp::Iterator<T>* it = iterable.iterator();  LOCAL (it);
        it->get (items);

        return items.size() == nbItems;
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_ */
