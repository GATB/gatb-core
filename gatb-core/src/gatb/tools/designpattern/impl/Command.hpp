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

/** \file Command.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of ICommand and ICommandDispatcher
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_

#include <gatb/tools/designpattern/api/ICommand.hpp>

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
    void dispatchCommands (std::vector<ICommand*>& commands, ICommand* post=0);

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
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_ */
