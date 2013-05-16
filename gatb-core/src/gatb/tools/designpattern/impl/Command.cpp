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

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/system/impl/System.hpp>

using namespace std;
using namespace gatb::core::tools::dp;

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/********************************************************************************/
class SynchronizerNull : public system::ISynchronizer
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
void SerialCommandDispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* post)
{
    for (std::vector<ICommand*>::iterator it = commands.begin(); it != commands.end(); it++)
    {
        if (*it != 0)  {  (*it)->use ();  (*it)->execute ();  (*it)->forget ();  }
    }

    /** We may have to do some post treatment. Note that we do it in the current thread. */
    if (post != 0)  {  post->use ();  post->execute ();  post->forget ();  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::ISynchronizer* SerialCommandDispatcher::newSynchro ()
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
ParallelCommandDispatcher::ParallelCommandDispatcher (size_t nbUnits) : _nbUnits(nbUnits)
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
void ParallelCommandDispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment)
{
    std::vector<system::IThread*> threads;
    for (std::vector<ICommand*>::iterator it = commands.begin(); it != commands.end(); it++)
    {
        threads.push_back (newThread(*it));
    }

    for (std::vector<system::IThread*>::iterator it = threads.begin(); it != threads.end(); it++)
    {
        (*it)->join ();
    }

    for (std::vector<system::IThread*>::iterator it = threads.begin(); it != threads.end(); it++)
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
system::ISynchronizer* ParallelCommandDispatcher::newSynchro ()
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
system::IThread* ParallelCommandDispatcher::newThread (ICommand* command)
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
void* ParallelCommandDispatcher::mainloop (void* data)
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
