/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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

/** \brief TO BE DONE...
 */
template <class Container, typename Type> class STLIterator : public Iterator<Type>
{
public:
    /** Constructor.
     * \param[in]  l : the list to be iterated
     */
    STLIterator (const Container& l)  : _l(l),_isDone(true) {}

    /** Destructor (here because of virtual methods). */
    virtual ~STLIterator ()  {}

    /** */
    void first()
    {
        _iter = _l.begin();
        _isDone = _iter == _l.end ();
        if (!_isDone)  { * this->_item = *_iter; }
    }

    /** */
    void next()
    {
        _iter++;
        _isDone = _iter == _l.end ();
        if (!_isDone)  { * this->_item = *_iter; }
    }

    /** */
    bool isDone() { return _isDone; }

    /** */
    Type& item()  { return * this->_item; }

private:

    /** List to be iterated. */
    Container _l;

    /** STL iterator we are going to wrap. */
    typename Container::iterator _iter;

    bool _isDone;
};

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
template <class Type> class ListIterator : public STLIterator<std::list<Type>, Type>
{
public:
    ListIterator (const std::list<Type>& l)  :  STLIterator<std::list<Type>, Type> (l)  {}
};

/********************************************************************************/

template <class Type> class VecIterator : public STLIterator<std::vector<Type>, Type>
{
public:
    VecIterator (const std::vector<Type>& l)  :  STLIterator<std::vector<Type>, Type> (l)  {}
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_WRAPPERS_HPP_ */
