/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Command.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of ICommand and ICommandDispatcher
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/misc/api/Macros.hpp>

#include <list>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** \brief Launches commands in current thread
 *
 * A dispatcher that uses the calling thread, so no parallelization.
 *
 * This implementation can be useful to process ICommand instances in an serial way
 * when it is required, while keeping an uniform API (ie. call dispatchCommands)
 * for running ICommand instances.
 */
class SerialCommandDispatcher : public ICommandDispatcher
{
public:

    /** Destructor (defined because of presence of virtual methods). */
    virtual ~SerialCommandDispatcher() {}

    /** \copydoc ICommandDispatcher::dispatchCommands */
    void dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment=0);

    /** \copydoc ICommandDispatcher::getExecutionUnitsNumber */
    size_t getExecutionUnitsNumber () { return 1; }

private:

    /** */
    system::ISynchronizer* newSynchro ();
};

/********************************************************************************/

/** \brief Launches commands in different threads
 *
 *  This implementation launches commands through different ICommandInvoker => parallelization.
 *  A provided ICommandInvokerFactory is used for creating ICommandInvoker instances.
 *
 *  This implementation of ICommandDispatcher is central in the PLAST design because it allows
 *  to uses all available cores. If one knows the number N of available cores on the computer,
 *  one has just to split some job by creating N ICommand instances and then just dispatch these
 *  commands through a ParallelCommandDispatcher: each command will be launched in a separated
 *  thread, and, thanks to the operating system architecture, each thread should be processed
 *  on an available core.
 *
 *  Note: it wouldn't be reasonable to use more ICommand instances than available cores.
 *  By default, if the number of dispatching units is not provided in the constructor of
 *  ParallelCommandDispatcher, it retrieves the number of available cores through the
 *  DefaultFactory::thread().getNbCores() method, and uses it as default value. This means
 *  that default constructor will use by default the whole CPU multicore power.
 */
class ParallelCommandDispatcher : public ICommandDispatcher
{
public:

    /** Constructor.
     * \param[in] nbUnits : number of threads to be used. If 0 is provided, one tries to guess the number of available cores.
     */
    ParallelCommandDispatcher (size_t nbUnits=0);

    /** \copydoc ICommandDispatcher::dispatchCommands */
    void dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment=0);

    /** \copydoc ICommandDispatcher::getExecutionUnitsNumber */
    size_t getExecutionUnitsNumber () { return _nbUnits; }

private:

    /** */
    system::ISynchronizer* newSynchro ();

    /** */
    system::IThread* newThread (ICommand* command);

    /** */
    static void* mainloop (void* data);

    /** Number of execution units to be used for command dispatching. */
    size_t _nbUnits;
};

/********************************************************************************/

class IteratorFunctor
{
public:

    IteratorFunctor ()  {}

    ~IteratorFunctor ()
    {
        /** We delete only synchronizers that have been deleted by the current instance. */
        for (std::list<Item>::iterator it = _synchros.begin(); it != _synchros.end(); it++)
        {
            if (this == it->first &&  it->second != 0)  { delete it->second; }
        }
    }

    system::ISynchronizer* newSynchro ()
    {
        system::ISynchronizer* synchro = system::impl::System::thread().newSynchronizer ();
        _synchros.push_back (Item(this,synchro));
        return synchro;
    }

private:

    typedef std::pair <IteratorFunctor*, system::ISynchronizer*> Item;

    std::list <Item>_synchros;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_ */
