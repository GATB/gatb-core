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

/** \file Iterator.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Smart Pointer Design Pattern definition
 *
 *  This file holds interfaces related to the Design Pattern Observer.
 */

#ifndef _GATB_CORE_DP_ITERATOR_HPP_
#define _GATB_CORE_DP_ITERATOR_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/system/api/ISmartPointer.hpp>

#include <vector>
#include <iostream>

/********************************************************************************/
namespace gatb  {
namespace core  {
/** \brief Tools package */
namespace tools {
/** \brief Design Patterns concepts */
namespace dp    {
/********************************************************************************/

/** \brief  Definition of the Design Pattern Iterator interface.
 *
 *    The Iterator concept is here reified as a template class that knows how to iterate some set of objects.
 *
 *  Actually, the interface has two ways for iterating instances:
 *    1- the 'classic' one in terms of Iterator Design Pattern (see first/next/isDone methods)
 *    2- a callback way (see 'iterate' methods) where some client provides a callback that will be called for each object
 *
 *  There may be good reasons for using one way or another. For instance, the first one may be easier to use by clients
 *  (no extra method to be defined) but may be less efficient because more methods calls will be carried out.
 *
 *  Note that 'iterate' has great subclassing potentiality; by default, its implementation relies on first/isDone/next/currentItem
 *  methods and so can't be faster that the 'classic' way. But specific implementations of Iterator may implement 'iterate'
 *  directly by using inner state of the class, and so without call to  first/isDone/next/currentItem methods.
 *
 *  Note that for optimization point of view, PLAST algorithm prefers to use the 'iterate' way with optimized implementations.
 *
 *  Sample of use:
 *  \code
 *   dp::Iterator<MyType*>* it = new MyIterator ();
 *   for (it->first(); ! it->isDone(); it->next() )
 *   {
 *      // retrieve the current item of some type
 *      MyType* item = it->currentItem ();
 *   }
 *  \endcode
 *
 *
 *  A remark on the STL: the Standard Template Library provides also the iterator concept. For instance, one could write:
 *  \code
 *   std::list<MyType*> l;
 *
 *   // we suppose here that we add some items to the list.
 *
 *   for (std::list<MyType*>::iterator it = l.begin(); it != l.end(); it++)
 *   {
 *      // retrieve the current item of some type
 *      MyType* item = *it;
 *   }
 *  \endcode
 *
 *  We can see the following differences with our own iterator.
 *   \li the iteration loop with the STL needs to know the iterated list, through the beginning iterator 'l.begin()'
 *   and the ending iterator 'l.end()'. In our case, we don't have any reference of the iterated container, we can see
 *   only a reference on the iterator itself.
 *   \li the iteration loop with the STL makes apparent the actual type of the container ('list' in the example), which
 *   is not the case for our iterator, where we can see only the type of the iterated items.
 *   \li the STL iterators can only iterate contents of containers, ie. some items already in memory. Our iterator concept
 *   is more general because it can iterate items that it builds on the fly, without having the whole items mounted in memory.
 *   A typical use is the parsing of a huge FASTA database (say 8 GBytes for instance). With a STL iterator, we should have
 *   the whole sequences lying in some containers before iterating them; with our iterator, it is enough to know a sequence
 *   at a given time. Note however that clients of such iterators must not reference a sequence that is transient by nature.
 *   In such a case, clients have to copy, if needed, the iterated sequence for some local work.
 *
 *  In brief, our iterator encapsulates more information that the STL iterator, allows to iterate potential containers (versus
 *  actual containers) and may be easier to use.
 *
 *  Moreover, we can use our iterator as a basis for other ways for iteration.
 */
template <class Item> class Iterator : public system::SmartPointer
{
public:

    /** */
    Iterator () : _item(&_default), _isRunning(false) {}

    /** Method that initializes the iteration. */
    virtual void first() = 0;

    /** Method that goes to the next item in the iteration.
     * \return status of the iteration
     */
    virtual void next()  = 0;

    /** Method telling whether the iteration is finished or not.
     * \return true if iteration is finished, false otherwise.
     */
    virtual bool isDone() = 0;

    /** Method that returns the current iterated item. Note that the returned type is the template type.
        \return the current item in the iteration.
    */
    virtual Item& item () = 0;

    /** Operator overload as a shortcut for the 'item' method. */
    Item* operator-> ()  { return &item(); }

    /** Operator overload as a shortcut for the 'item' method. */
    Item& operator*  ()  { return item();  }

    /** Another way to iterate: push model, ie a functor is called for each item. */
    template <typename Functor> void iterate (const Functor& f)   {  for (first(); !isDone(); next())  { f (item()); }  }

    /** Get a reference on the object to be configured as the currently iterated item.
     * \param[in] i : object to be referred. */
    virtual void setItem (Item& i)  {  _item = &i;  }

    /** Retrieve some iterated items in a vector.
     * \param[in] current : vector to be filled with iterated items. May be resized if not enough items available
     * \return true if the iteration is not finished, false otherwise. */
    bool get (std::vector<Item>& current)
    {
        size_t i=0;
        for (i=0; i<current.size(); i++)
        {
            setItem (current[i]);

            if (_isRunning == false)  { first ();  _isRunning=true; }
            else                      { next  ();                   }

            if (isDone())
            {
                current.resize (i);
                return false;
            }
        }
        return true;
    }

    /** Reset the iterator. */
    virtual void reset ()
    {
        _isRunning = false;
        _item      = &_default;
    }

protected:
    Item* _item;

private:
    Item  _default;
    bool  _isRunning;
};

/********************************************************************************/

template<typename T>
class ISmartIterator : public Iterator<T>
{
public:
    virtual ~ISmartIterator() {}
    virtual u_int64_t size () const = 0;
    virtual u_int64_t rank () const = 0;
};

/********************************************************************************/

/** \brief Interface for listening to iteration progress.
 *
 * This interface is intended to be notified by some progress job, and in particular
 * to the progression of an iteration with an Iterator instance.
 *
 * It defines 3 methods:
 *      - init   : method called just before the beginning of the iteration
 *      - finish : method called just after the end of the iteration
 *      - inc    : method called during the iteration; the provided arguments
 *                 represents the number of iterations before the previous call
 *
 * Actually, this is a little more than just an interface since it provides empty
 * implementations for the 3 methods; this will ease the development to clients who
 * wants only to get 'inc' notifications but are not interested to do specific actions
 * at the beginning and the end of the iteration.
 *
 * \see SubjectIterator
 */
class IteratorListener : public system::SmartPointer
{
public:

    /** Destructor. */
    virtual ~IteratorListener ()  {}

    /** Initialization of the object. */
    virtual void init () {}

    /** Finish the progress information. */
    virtual void finish () {}

    /** Increase the number of currently done tasks. */
    virtual void inc (u_int64_t ntasks_done) {}

    /** Associate a message to the listener.
     * \param[in] format : printf-like argument for building the message */
    virtual void setMessage (const char* format, ...)  {}

    /** Set the current number of tasks done.
     * \param[in] ntasks_done :  sets the current number of job done. */
    virtual void set (u_int64_t ntasks_done) {}
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_HPP_ */
