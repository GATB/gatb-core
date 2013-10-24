/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
ParallelDispatcher::ParallelDispatcher (size_t nbUnits) : _nbUnits(nbUnits)
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
size_t ParallelDispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment)
{
    TIME_START (ti, "compute");

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
system::ISynchronizer* ParallelDispatcher::newSynchro ()
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
system::IThread* ParallelDispatcher::newThread (ICommand* command)
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
void* ParallelDispatcher::mainloop (void* data)
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
