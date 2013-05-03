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

#ifdef __APPLE__

#include <gatb/system/impl/ThreadMacos.hpp>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <unistd.h>

using namespace std;

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
class ThreadMacos : public IThread
{
public:
	ThreadMacos (void* (mainloop) (void*), void* data)  { pthread_create (&_thread, NULL,  mainloop, data); }
    ~ThreadMacos ()  { /* pthread_detach (_thread); */  }
    void join ()     { pthread_join   (_thread, NULL);  }
private:
    pthread_t  _thread;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
class SynchronizerMacos : public ISynchronizer
{
public:
	SynchronizerMacos ()            {  pthread_mutex_init (&_mutex, NULL);  }
    virtual ~SynchronizerMacos()    {  pthread_mutex_destroy (&_mutex);     }

    void   lock ()  { pthread_mutex_lock   (&_mutex); }
    void unlock ()  { pthread_mutex_unlock (&_mutex); }

private:
    pthread_mutex_t  _mutex;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThread* ThreadFactoryMacos::newThread (void* (*mainloop) (void*), void* data)
{
    return new ThreadMacos (mainloop, data);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ISynchronizer* ThreadFactoryMacos::newSynchronizer (void)
{
    return new SynchronizerMacos ();
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* __LINUX__ */
