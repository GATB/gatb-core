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

/** \file IteratorWrappers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Some wrappers classes for our Iterator design
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_WRAPPERS_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_WRAPPERS_HPP_

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <list>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** \brief Iterator that loops over std::list
 *
 *  This class is a wrapper between the STL list implementation and our own Iterator
 *  abstraction.
 *
 *  Note that the class is still a template one since we can iterate on list of anything.
 *
 *  \code
 *  void foo ()
 *  {
 *      list<int> l;
 *      l.push_back (1);
 *      l.push_back (2);
 *      l.push_back (4);
 *
 *      ListIterator<int> it (l);
 *      for (it.first(); !it.isDone(); it.next())       {   cout << it.currentItem() << endl;   }
 *  }
 *  \endcode
 */
template <class T1> class ListIterator : public Iterator<T1>
{
public:
    /** Constructor.
     * \param[in]  l : the list to be iterated
     */
    ListIterator (const std::list<T1>& l)  : _l(l) {}

    /** Destructor (here because of virtual methods). */
    virtual ~ListIterator ()  {}

    /** \copydoc Iterator<T1>::first */
    void first()  {  _iter = _l.begin(); }

    /** \copydoc Iterator<T1>::next */
    void next()   { _iter++; }

    /** \copydoc Iterator<T1>::isDone */
    bool isDone() { return _iter == _l.end (); }

    /** \copydoc Iterator<T1>::item */
    T1& item()  { return *_iter; }

private:

    /** List to be iterated. */
    std::list<T1> _l;

    /** STL iterator we are going to wrap. */
    typename std::list<T1>::iterator _iter;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_WRAPPERS_HPP_ */
