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

/** \file ThreadLinux.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Thread management for Linux
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_LINUX_THREAD_HPP_
#define _GATB_CORE_SYSTEM_IMPL_LINUX_THREAD_HPP_

#include <gatb/system/api/IThread.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
namespace impl      {
/********************************************************************************/

/** \brief Factory that creates IThread instances.
 *
 *  Thread creation needs merely the main loop function that will be called.
 */
class ThreadFactoryLinux : public IThreadFactory
{
public:
    /** \copydoc IThreadFactory::newThread */
    IThread* newThread (void* (*mainloop) (void*), void* data);

    /** \copydoc IThreadFactory::newSynchronizer */
    ISynchronizer* newSynchronizer (void);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IMPL_LINUX_THREAD_HPP_ */
