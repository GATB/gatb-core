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

#ifndef _GATB_CORE_SYSTEM_ITHREAD_HPP_
#define _GATB_CORE_SYSTEM_ITHREAD_HPP_

#include <gatb/system/api/types.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/********************************************************************************/

/** \brief Define what a thread is.
 *
 * Definition of a thread in an OS independent fashion. See below how to create a thread through a factory.
 */
class IThread
{
public:

    /** Wait the end of the thread. */
    virtual void join () = 0;

    /** Destructor. */
    virtual ~IThread () {}
};

/********************************************************************************/

/** \brief Define a synchronization abstraction
 *
 *  This is an abstraction layer of what we need for handling synchronization.
 *  Actual implementations may use mutex for instance.
 */
class ISynchronizer
{
public:

    /** Lock the synchronizer. */
    virtual void   lock () = 0;

    /** Unlock the synchronizer. */
    virtual void unlock () = 0;

    /** Destructor. */
    virtual ~ISynchronizer () {}
};

/********************************************************************************/

/** \brief Factory that creates IThread instances.
 *
 *  Thread creation needs merely the main loop function that will be called.
 *
 *  Note the method that can return the number of cores in case a multi-cores
 *  architecture is used. This is useful for automatically configure PLAST for
 *  using the maximum number of available cores for speeding up the algorithm.
 */
class IThreadFactory
{
public:

    /** Creates a new thread.
     * \param[in] mainloop : the function the thread shall execute
     * \param[in] data :  data provided to the mainloop when launched
     * \return the created thread.
     */
    virtual IThread* newThread (void* (*mainloop) (void*), void* data) = 0;

    /** Creates a new synchronization object.
     * \return the created ISynchronizer instance
     */
    virtual ISynchronizer* newSynchronizer (void) = 0;

    /** Destructor. */
    virtual ~IThreadFactory ()  {}
};

/********************************************************************************/

/** \brief Tool for locally managing synchronization.
 *
 *  Instances of this class reference a ISynchronizer instance. When created, they lock
 *  their referred ISynchronizer and when destroyed, they unlock it.
 *
 *  For instance, it is a convenient way to lock/unlock a ISynchronizer instance in the scope
 *  of a statements block (where the LocalSynchronizer instance lives), a particular case being
 *  the statements block of a method.
 *
 *  Code sample:
 *  \code
 *  void sample (ISynchronizer* synchronizer)
 *  {
 *      // we create a local synchronizer from the provided argument
 *      LocalSynchronizer localsynchro (synchronizer);
 *
 *      // now, all the statements block is locked for the provided synchronizer
 *      // in other word, this sample function is protected against concurrent accesses.
 *  }
 *  \endcode
 */
class LocalSynchronizer
{
public:

    /** Constructor.
     * \param[in] ref : the ISynchronizer instance to be controlled.
     */
    LocalSynchronizer (ISynchronizer* ref) : _ref(ref)  {  if (_ref)  { _ref->lock (); }  }

    /** Destructor. */
    ~LocalSynchronizer ()  {  if (_ref)  {  _ref->unlock (); } }

private:

    /** The referred ISynchronizer instance. */
    ISynchronizer* _ref;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_ITHREAD_HPP_ */
