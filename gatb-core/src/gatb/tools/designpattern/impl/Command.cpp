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

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

using namespace std;
using namespace gatb::core::tools::dp;

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** */
class CommandStartSynchro : public ICommand, public system::SmartPointer
{
public:

    CommandStartSynchro (ICommand* ref, system::ISynchronizer* synchro) : _ref(0), _synchro(synchro)  { setRef(ref); }
    ~CommandStartSynchro ()  { setRef(0); }

    void execute ()
    {
        if (_ref && _synchro)
        {
            /** We lock/unlock the synchronizer. */
            _synchro->lock();  _synchro->unlock ();

            /** We execute the delegate command. */
            _ref->execute();
        }
    }

private:
    ICommand* _ref;
    void setRef (ICommand* ref)  { SP_SETATTR(ref); }

    system::ISynchronizer* _synchro;
};

/********************************************************************************/
class SynchronizerNull : public system::ISynchronizer, public system::SmartPointer
{
public:
    void   lock ()  {}
    void unlock ()  {}
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t SerialDispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* post)
{
    TIME_START (ti, "compute");

    for (std::vector<ICommand*>::iterator it = commands.begin(); it != commands.end(); it++)
    {
        if (*it != 0)  {  (*it)->use ();  (*it)->execute ();  (*it)->forget ();  }
    }

    /** We may have to do some post treatment. Note that we do it in the current thread. */
    if (post != 0)  {  post->use ();  post->execute ();  post->forget ();  }

    TIME_STOP (ti, "compute");

    return ti.getEntryByKey("compute");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::ISynchronizer* SerialDispatcher::newSynchro ()
{
    return new SynchronizerNull();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Dispatcher::Dispatcher (size_t nbUnits, size_t groupSize) : _nbUnits(nbUnits), _groupSize(groupSize)
{
    if (_nbUnits==0)  { _nbUnits = system::impl::System::info().getNbCores(); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t Dispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment)
{
    TIME_START (ti, "compute");

    system::IThreadGroup* threadGroup = system::impl::ThreadGroup::create ();

    for (std::vector<ICommand*>::iterator it = commands.begin(); it != commands.end(); it++)
    {
        /** We add the thread to the group. */
        threadGroup->add (mainloop, new CommandStartSynchro (*it, threadGroup->getSynchro()) );
    }

    /** We start the group. */
    threadGroup->start ();

    /** Some cleanup. */
    system::impl::ThreadGroup::destroy (threadGroup);

    TIME_STOP (ti, "compute");

    /** We return the result. */
    return ti.getEntryByKey("compute");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::ISynchronizer* Dispatcher::newSynchro ()
{
    return system::impl::System::thread().newSynchronizer();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::IThread* Dispatcher::newThread (ICommand* command)
{
    return system::impl::System::thread().newThread (mainloop, command);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void* Dispatcher::mainloop (void* data)
{
    ICommand* cmd = (ICommand*) data;

    if (cmd != 0)
    {
        cmd->use ();
        cmd->execute();
        cmd->forget ();
    }

    return 0;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
