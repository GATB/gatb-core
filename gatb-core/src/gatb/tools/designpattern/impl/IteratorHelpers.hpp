/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file IteratorHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Some helper classes for iteration
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_HELPERSHPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_HELPERSHPP_

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <set>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** \brief Null implementation of the Iterator interface.
 *
 * This implementation merely iterates over nothing. It may be useful to use such
 * an instance when we have to provide an Iterator instance but with nothing to iterate.
 */
template <class Item> class NullIterator : public Iterator<Item>
{
public:

    /** \copydoc  Iterator::first */
    void first() {}

    /** \copydoc  Iterator::next */
    void next()  {}

    /** \copydoc  Iterator::isDone */
    bool isDone() { return true; }

    /** \copydoc  Iterator::item */
    Item& item ()  { static Item item;  return item; }

    /** Destructor. */
    ~NullIterator() {}
};

/********************************************************************************/

/** \brief Iterator over two iterators.
 *
 * We define a "product" iterator for two iterators, i.e. it will loop each possible
 * couple of the two provided iterators.
 *
 *  It is useful for having only one loop instead of two loops. Note however that it
 *  may still be more efficient to have two loops. The CartesianIterator is just here
 *  for easing the product iteration on small sets.
 *
 *  \code
 *  // We declare a Bank instance.
 *  Bank b (filename);
 *
 *  // We declare two iterators on the bank
 *  Bank::Iterator it1 (b);
 *  Bank::Iterator it2 (b);
 *
 *  // We build a product iterator from the two iterators
 *  CartesianIterator<Sequence,Sequence> prodIt (it1, it2);
 *
 *  for (prodIt.first(); !prodIt.isDone(); prodIt.next())
 *  {
 *      printf ("'%s'\n'%s'\n\n", prodIt->first.getData(), prodIt->second.getData());
 *  }
 *  \endcode
 */
template <class T1, class T2> class CartesianIterator : public Iterator < std::pair<T1,T2> >
{
public:

    /** Constructor.
     * \param[in] it1 : first iterator.
     * \param[in] it2 : second iterator.
     */
    CartesianIterator (Iterator<T1>& it1, Iterator<T2>& it2)  : _it1(it1), _it2(it2), _isDone(false)  {  first(); }

    /** Destructor. */
    virtual ~CartesianIterator ()  {}

    /** \copydoc Iterator::first */
    void first()
    {
        /** We go to the beginning of the two iterators. */
        _it1.first();
        _it2.first();

        /** We use a specific attribute to keep track of the finish status of the iterator.
         *  This is merely an optimization in order not to call too often the "isDone" method
         *  on the two iterators.  */
        _isDone = false;
    }

    /** \copydoc Iterator::next */
    void next()
    {
        /** We go to the next item of the second iterator. */
        _it2.next ();

        /** We check that it is not done. */
        if (_it2.isDone())
        {
            /** We go to the next item of the first iterator. */
            _it1.next ();

            if (! _it1.isDone() )
            {
                /** The first iterator is not done, we can reset the second to its beginning. */
                _it2.first ();
            }
            else
            {
                /** The first iterator is also done, the product iterator is therefore done. */
                _isDone = true;
            }
        }
    }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _isDone; }

    /** \copydoc Iterator::item */
    std::pair<T1,T2>& item ()
    {
        _current.first  = _it1.item();
        _current.second = _it2.item();

        return _current;
    }

private:

    /** First iterator. */
    Iterator<T1>& _it1;

    /** Second iterator. */
    Iterator<T2>& _it2;

    /** Current item in the iteration. */
    std::pair<T1,T2> _current;

    /** Tells whether the iteration is finished or not. */
    bool _isDone;
};

/********************************************************************************/

/** \brief Iterator over two iterators.
 *
 * We define a an iterator for two iterators, by iterating each of the two iterators
 * and providing a pair of the two currently iterated items.
 */
template <class T1, class T2> class PairedIterator : public Iterator < std::pair<T1,T2> >
{
public:

    /** Constructor.
     * \param[in] it1 : first iterator.
     * \param[in] it2 : second iterator.
     */
    PairedIterator (Iterator<T1>& it1, Iterator<T2>& it2)  : _it1(it1), _it2(it2) { }

    /** Destructor. */
    virtual ~PairedIterator ()  {}

    /** \copydoc Iterator::first */
    void first()  {  _it1.first();  _it2.first();  }

    /** \copydoc Iterator::next */
    void next()  {  _it1.next ();  _it2.next ();  }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _it1.isDone() || _it2.isDone(); }

    /** \copydoc Iterator::item */
    std::pair<T1,T2>& item ()
    {
        _current.first  = _it1.item();
        _current.second = _it2.item();

        return _current;
    }

private:

    /** First iterator. */
    Iterator<T1>& _it1;

    /** Second iterator. */
    Iterator<T2>& _it2;

    /** Current item in the iteration. */
    std::pair<T1,T2> _current;
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
class IteratorListener : public SmartPointer
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

    /** Associate a message to the listener. */
    virtual void setMessage (const char* format, ...)  {}
};

/********************************************************************************/

/** \brief Factorization of code for the subject part of the Observer pattern.
 */
class AbstractSubjectIterator
{
public:

    /** Constructor. */
    AbstractSubjectIterator () : _hasListeners(false), _isStarted(false) {}

    /** Destructor. */
    ~AbstractSubjectIterator ()
    {
        /** We remove all observers. */
        for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
        {
            (*it)->forget();
        }
    }

    /** Add an observer to the iterator. Such an observer is provided as a functor.
     * \param[in] f : functor to be subscribed to the iterator notifications.
     */
    void addObserver    (IteratorListener* f)
    {
        if (f != 0)
        {
            f->use ();
            _listeners.insert (f);
            _hasListeners=true;
        }
    }

    /** Remove an observer from the iterator. Such an observer is provided as a functor.
     * \param[in] f : functor to be unsubscribed from the iterator notifications.
     */
    void removeObserver (IteratorListener& f)
    {
        /** We look whether the given functor is already known. */
        std::set<IteratorListener*>::iterator lookup = _listeners.find (&f);
        if (lookup != _listeners.end())
        {
            (*lookup)->forget();
            _listeners.erase (lookup);
            _hasListeners = _listeners.empty() == false;
        }
    }

protected:

    /** Notify all the subscribed functors.
     * \param[in] current : number of currently iterated items during the iteration. */
    void notifyInc (u_int64_t current)
    {
        if (_isStarted == true)
        {
            /** We call each subscribing functor. */
            for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
            {
                /** Not very pretty syntax, but that works. */
               (*it)->inc (current);
            }
        }
    }

    /** Notify all the subscribed functors about the start of the iteration. */
    void notifyInit ()
    {
        if (_isStarted == false)
        {
            _isStarted = true;

            /** We call each subscribing functor. */
            for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
            {
                (*it)->init ();
            }
        }
    }

    /** Notify all the subscribed functors about the start of the iteration. */
    void notifyFinish ()
    {
        if (_isStarted == true)
        {
            _isStarted = false;

            /** We call each subscribing functor. */
            for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
            {
                (*it)->finish ();
            }
        }
    }

private:

    std::set<IteratorListener*> _listeners;
    bool                        _hasListeners;
    bool                        _isStarted;
};

/********************************************************************************/

/** \brief Iterator that notifies some listener during iteration.
 *
 * Implementation note: we have to keep reference (through pointers) on functors
 * because we want them to be notified and not a copy of them.
 *
 * Note also that we don't allow to have twice the same observer (we use a set as
 * observers container).
 *
 * \snippet iterators4.cpp  snippet1_SubjectIterator
 */
template <class Item> class SubjectIterator : public Iterator<Item>, public AbstractSubjectIterator
{
public:

    /** Constructor
     * \param[in] ref : the referred iterator
     * \param[in] modulo : notifies every 'modulo' time */
    SubjectIterator (Iterator<Item>* ref, u_int32_t modulo)  : _ref(0), _modulo(modulo==0 ? 1 : modulo), _current(0)
    {
        /** We set the reference. */
        setRef (ref);
    }

    /** Destructor. */
    ~SubjectIterator ()
    {
        /** We release the reference. */
        setRef (0);
    }

    /** \copydoc Iterator::first */
    void first ()
    {
        notifyInit ();
        _current = 0;
        _ref->first ();
    }

    /** \copydoc Iterator::isDone */
    bool isDone ()
    {
        bool res = _ref->isDone();  if (res)  { notifyFinish(); }  return res;
    }

    /** \copydoc Iterator::next */
    void next ()
    {
        _ref->next ();
        if ((_current % _modulo) == 0)  { notifyInc (_current);  _current=0; }
        _current++;
    }

    /** \copydoc Iterator::item */
    Item& item ()  { return _ref->item(); }

    /** */
    void setItem (Item& current)  { _ref->setItem(current); }

    /** */
    void reset ()  { _ref->reset(); }

private:

    Iterator<Item>* _ref;
    void setRef (Iterator<Item>* ref)  { SP_SETATTR(ref); }

    u_int64_t _current;
    u_int64_t _modulo;
};

/********************************************************************************/
/** \brief Iterator that gather two iterators in a single loop.
 *
 * This iterator is equivalent to two iterators with one outer loop and one inner loop.
 */
template <typename T1, typename T2, typename Update> class CompoundIterator : public Iterator<T2>
{
public:
    /** Constructor.
     * \param[in] it1 : iterator for the outer loop.
     * \param[in] it2 : iterator for the inner loop.
     * \param[in] update : functor for updating inner loop with new item from outer loop
     */
    CompoundIterator (Iterator<T1>&  it1, Iterator<T2>& it2, const Update& update)  : _it1(it1), _it2(it2), _update(update) {}

    /** \copydoc Iterator::first */
    void first ()
    {
        _it1.first();
        if (!_it1.isDone())
        {
            _update (&_it2, &_it1.item());
            _it2.first ();
        }
    }

    /** \copydoc Iterator::next */
    void next ()
    {
        _it2.next ();
        if (_it2.isDone())
        {
            _it1.next();
            if (!_it1.isDone())
            {
                _update (&_it2, &_it1.item());
                _it2.first ();
            }
        }
    }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _it1.isDone(); }

    /** \copydoc Iterator::item */
    T2& item()  { return _it2.item(); }

private:
    Iterator<T1>&  _it1;
    Iterator<T2>&  _it2;
    Update         _update;
};

/********************************************************************************/
/** \brief Iterator that truncate the number of iterations if needed.
 *
 * This iterator iterates a referred iterator and will finish:
 *      - when the referred iterator is over
 *   or - when a limit number of iterations is reached.
 */
template <class Item> class TruncateIterator : public Iterator<Item>
{
public:

    /** Constructor.
     * \param[in] ref : the referred iterator
     * \param[in] limit : the maximal number of iterations. */
    TruncateIterator (Iterator<Item>& ref, u_int64_t limit) : _ref(ref), _limit(limit), _currentIdx(0) {}

    /** \copydoc  Iterator::first */
    void first() { _currentIdx=0;  _ref.first(); }

    /** \copydoc  Iterator::next */
    void next()  { _currentIdx++;  _ref.next(); }

    /** \copydoc  Iterator::isDone */
    bool isDone() { return _ref.isDone() || _currentIdx >= _limit; }

    /** \copydoc  Iterator::item */
    Item& item ()  { return _ref.item(); }

    /** \copydoc  Iterator::setItem */
    void setItem (Item& i)  { _ref.setItem(i); }

private:

    Iterator<Item>& _ref;
    u_int64_t       _limit;
    u_int64_t       _currentIdx;
};

/********************************************************************************/
/** \brief Iterator that filters out some iterated items
 *
 * This iterator iterates a referred iterator and will filter out some items according
 * to a functor provided at construction.
 */
template <class Item, typename Filter> class FilterIterator : public Iterator<Item>
{
public:

    /** Constructor.
     * \param[in] ref : the referred iterator
     * \param[in] filter : the filter on items. Returns true if item is kept, false otherwise. */
    FilterIterator (Iterator<Item>& ref, Filter& filter) : _ref(ref), _filter(filter)  {}

    /** \copydoc  Iterator::first */
    void first() { _ref.first(); while (!isDone() && _filter(item())==false) { _ref.next(); }  }

    /** \copydoc  Iterator::next */
    void next()  { _ref.next();  while (!isDone() && _filter(item())==false) { _ref.next(); }  }

    /** \copydoc  Iterator::isDone */
    bool isDone() { return _ref.isDone();  }

    /** \copydoc  Iterator::item */
    Item& item ()  { return _ref.item(); }

    /** \copydoc  Iterator::setItem */
    void setItem (Item& i)  { _ref.setItem(i); }

private:

    Iterator<Item>& _ref;
    Filter&         _filter;
};

/********************************************************************************/

template <class Item> class VectorIterator : public Iterator<Item>
{
public:
    /** */
    VectorIterator () : _idx(0), _nb (0)  {}

    /** */
    virtual ~VectorIterator () {}

    /** \copydoc  Iterator::first */
    void first()  {  _idx = -1;  next ();  }

    /** \copydoc  Iterator::next */
    void next()  { ++_idx;  if (_idx < _nb ) { this->_item = &(_items[_idx]); }  }

    /** \copydoc  Iterator::isDone */
    bool isDone() {   return _idx >= _nb;  }

    /** \copydoc  Iterator::item */
    Item& item ()  { return *(this->_item); }

protected:
    std::vector<Item> _items;
    int32_t           _idx;
    int32_t           _nb;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_HELPERSHPP_ */
