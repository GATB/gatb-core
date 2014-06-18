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

/** \file System.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Entry point class for accessing operating system operations
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_SYSTEM_HPP_
#define _GATB_CORE_SYSTEM_IMPL_SYSTEM_HPP_

/********************************************************************************/

#include <gatb/system/impl/MemoryCommon.hpp>
#include <gatb/system/impl/TimeCommon.hpp>
#include <gatb/system/impl/SystemInfoCommon.hpp>
#include <gatb/system/impl/FileSystemLinux.hpp>
#include <gatb/system/impl/FileSystemMacos.hpp>
#include <gatb/system/impl/ThreadLinux.hpp>
#include <gatb/system/impl/ThreadMacos.hpp>
#include <map>
#include <functional>
#include <list>
#include <vector>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {

/********************************************************************************/

#if 0
    #define MALLOC      System::memory().malloc
    #define CALLOC      System::memory().calloc
    #define REALLOC     System::memory().realloc
    #define FREE        System::memory().free
#else
    #define MALLOC      malloc
    #define CALLOC      calloc
    #define REALLOC     realloc
    #define FREE        free
#endif

/********************************************************************************/

 /** \brief Entry point providing access to operating system resources.
  *
  *  The IResource class provides a unique entry point for accessing different kinds
  *  of operating system resources (threads, time, file, etc...).
  *
  *  Normally, the client should use such an instance for getting OS resources
  *  instead of directly call system defendant functions. This is important because it
  *  will ease the build of client tools for different OS/architecture; the only thing
  *  to do is to use a specific instance of IResource for matching the correct OS.
  */
class System
 {
 public:

    /********************************************************************************/
    /** Access for info methods. */
    static ISystemInfo&  info  ()
    {
#ifdef __linux__
        static SystemInfoLinux instance;  return instance;
#endif

#ifdef __APPLE__
        static SystemInfoMacos instance;  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
    }

    /********************************************************************************/
     /** Access for time methods. */
     static ITime&  time ()
     {
#ifdef __linux__
        static TimeSystem instance (ITime::MSEC);  return instance;
#endif

#ifdef __APPLE__
        static TimeSystem instance (ITime::MSEC);  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

     /********************************************************************************/
     /** Access for file methods. */
     static IFileSystem&     file    ()
     {
#ifdef __linux__
        static FileSystemLinux instance;  return instance;
#endif

#ifdef __APPLE__
        static FileSystemMacos instance;  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

     /********************************************************************************/
     /** Access for memory methods. */
     static IMemory&         memory  ()
     {
#ifdef __linux__
        static MemoryCommon instance (MemoryAllocatorStdlib::singleton(), MemoryOperationsCommon::singleton());  return instance;
#endif

#ifdef __APPLE__
        static MemoryCommon instance (MemoryAllocatorStdlib::singleton(), MemoryOperationsCommon::singleton());  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

     /********************************************************************************/
     /** Access for thread methods. */
     static IThreadFactory&  thread  ()
     {
#ifdef __linux__
        static ThreadFactoryLinux instance;  return instance;
#endif

#ifdef __APPLE__
        static ThreadFactoryMacos instance;  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }
};

/********************************************************************************/

class ThreadGroup : public IThreadGroup, public system::SmartPointer
{
public:

    static IThreadGroup* create ();
    static void destroy (IThreadGroup* thr);

    static IThreadGroup* find (IThread::Id id);

    void add (void* (*mainloop) (void*), void* data);

    void start ();

    ISynchronizer* getSynchro()  { return _startSynchro; }

    size_t size() const { return _threads.size(); }

    IThread* operator[] (size_t idx)  { return _threads[idx]; }

    /** */
    void addException (system::Exception e)
    {
        LocalSynchronizer synchro (_startSynchro);
        _exceptions.push_back (e);
    }
    bool hasExceptions() const { return _exceptions.empty() == false; }

    Exception getException () const  { return ExceptionComposite(_exceptions); }

private:

    ThreadGroup ();
    ~ThreadGroup();

    std::vector<IThread*>  _threads;
    system::ISynchronizer* _startSynchro;

    static std::list<ThreadGroup*> _groups;

    std::list<system::Exception> _exceptions;
};

/********************************************************************************/
template<typename T> class ThreadObject
{
public:

    ThreadObject (const T& object = T()) : _object(object), _isInit(false), _synchro(0)
    {
        _synchro = system::impl::System::thread().newSynchronizer();
    }

    ~ThreadObject()
    {
        if (_synchro)  { delete _synchro; }

        for (typename std::map<IThread::Id,T*>::iterator it = _map.begin(); it != _map.end(); it++)  { delete it->second; }
    }

    T& operator () ()
    {

        if (_isInit == false)
        {
            LocalSynchronizer ls (_synchro);
            if (_isInit == false)
            {
                /** We look for the TreadGroup if any. */
                IThreadGroup* group = ThreadGroup::find (System::thread().getThreadSelf());

                if (group)
                {
                    for (size_t i=0; i<group->size(); i++)
                    {
                        IThread* thread = (*group)[i];
                        T* newObject = new T(this->_object);
                        this->_map[thread->getId()] = newObject;
                        this->_vec.push_back(newObject);
                    }
                }
                else
                {
                    std::cout << "ThreadObject::operator()  CAN'T HAPPEN......" << std::endl;
                }
                _isInit = true;
            }
        }

        return *(_map[System::thread().getThreadSelf()]);
    }

    void terminate ()  { }

    /** */
    template<typename Functor> void foreach (const Functor& fct)
    {  for (typename std::map<IThread::Id,T*>::iterator it = _map.begin(); it != _map.end(); it++)  { fct (*(it->second)); } }

    T* operator-> ()  { return &_object; }

    T& operator* ()  { return _object; }

    /** */
    size_t size() const { return _vec.size(); }

    /** */
    T& operator[] (size_t idx) { return *(_vec[idx]); }

private:
    std::map<IThread::Id,T*> _map;
    T _object;

    std::vector<T*> _vec;

    bool _isInit;
    system::ISynchronizer* _synchro;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_SYSTEM_IMPL_SYSTEM_HPP_ */
