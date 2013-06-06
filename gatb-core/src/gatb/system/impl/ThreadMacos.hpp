/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file ThreadMacos.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Thread management for Macos
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_MACOS_THREAD_HPP_
#define _GATB_CORE_SYSTEM_IMPL_MACOS_THREAD_HPP_

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
class ThreadFactoryMacos : public IThreadFactory
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

#endif /* _GATB_CORE_SYSTEM_IMPL_MACOS_THREAD_HPP_ */
