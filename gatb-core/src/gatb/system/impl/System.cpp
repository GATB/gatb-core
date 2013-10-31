/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <algorithm>

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

std::list<ThreadGroup*> ThreadGroup::_groups;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ThreadGroup::ThreadGroup()
:  _startSynchro(0)
{
    _startSynchro = System::thread().newSynchronizer();
    if (_startSynchro)  { _startSynchro->lock(); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ThreadGroup::~ThreadGroup()
{
    if (_startSynchro)  { delete _startSynchro; }

    /** We delete each thread. */
    for (std::vector<system::IThread*>::iterator it = _threads.begin(); it != _threads.end(); it++)
    {
        delete (*it);
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::add (void* (*mainloop) (void*), void* data)
{
    IThread* thr = System::thread().newThread (mainloop, data);

    _threads.push_back (thr);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::start ()
{
    /** We unlock the synchronizer => all threads of the group begin at the same time. */
    if (_startSynchro)  { _startSynchro->unlock(); }

    /** We join each thread. */
    for (std::vector<system::IThread*>::iterator it = _threads.begin(); it != _threads.end(); it++)
    {
        (*it)->join ();
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThreadGroup* ThreadGroup::create ()
{
    ThreadGroup* tg = new ThreadGroup;
    _groups.push_back(tg);
    return tg;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::destroy (IThreadGroup* thr)
{
    /** We look for the provided thread group. */
    for (std::list<ThreadGroup*>::iterator it = _groups.begin(); it != _groups.end(); it++)
    {
        if (*it == thr)
        {
            _groups.erase (it);
            delete thr;
            break;
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThreadGroup* ThreadGroup::find (IThread::Id id)
{
    /** We look for the provided thread group. */
    for (std::list<ThreadGroup*>::iterator it = _groups.begin(); it != _groups.end(); it++)
    {
        ThreadGroup* group = *it;

        for (std::vector<IThread*>::iterator itThread = group->_threads.begin();  itThread != group->_threads.end(); itThread++)
        {
            if ((*itThread)->getId() == id)  { return group; }
        }
    }
    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::foreach (const std::function<void (IThread*)>& fct)
{
    std::for_each (_threads.begin(), _threads.end(), fct);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
