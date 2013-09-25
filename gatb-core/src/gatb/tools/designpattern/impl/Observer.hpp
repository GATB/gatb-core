/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Observer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of INode interface
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_OBSERVER_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_OBSERVER_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/IObserver.hpp>
#include <list>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/


/** \brief Class that notifies potential observers.
 *
 * The main purpose of this class is to manage the set of IObservers attached to the subject.
 *
 * Then, classes that want subject-like behavior can inherit from Subject or have a Subject
 * attribute.
 *
 * \see ISubject
 */
class Subject : public ISubject
{
public:

    /** Constructor. */
    Subject ();

    /** Constructor.
     * \param[in] interface : the identifier of the subject. */
    Subject (const InterfaceId& interface);

    /** Destructor. */
    virtual ~Subject();

    /** Returns the identifier of the subject. */
    InterfaceId getInterface ()  { return _interface; }

    /** Attach an observer to this subject.
     * \param[in] observer : the observer to be attached.
     */
    void addObserver    (IObserver* observer);

    /** Detach an observer from this subject.
     * \param[in] observer : the observer to be detached.
     */
    void removeObserver (IObserver* observer);

    /** Notify observers by sending a EventInfo instance.
     * \param[in]  event : the information to be sent to the observers.
     */
    void notify (EventInfo* event);

private:

    /** Identifier of the subject. */
    InterfaceId             _interface;

    /** Optimization attribute. */
    IObserver*              _singleObserver;

    /** List of observers attached to the subject. */
    std::list<IObserver*>*  _observers;
};


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_CELL_HPP_ */
