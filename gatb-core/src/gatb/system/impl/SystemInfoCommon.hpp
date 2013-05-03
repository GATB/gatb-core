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

/** \file SystemInfoCommon.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementations common to various OS.
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_SYSTEM_COMMON_HPP_
#define _GATB_CORE_SYSTEM_IMPL_SYSTEM_COMMON_HPP_

/********************************************************************************/

#include <gatb/system/api/ISystemInfo.hpp>
#include <gatb/system/api/Exception.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {
/********************************************************************************/

class SystemInfoCommon : public ISystemInfo
{
public:

    /** \copydoc ISystemInfo::version */
    std::string version ()  { return std::string(__DATE__); }
};

/********************************************************************************/

class SystemInfoLinux : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores ();

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName ();

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal ();

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed ();

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers ();
};

/********************************************************************************/

class SystemInfoMacos : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores ();

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName ();

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal ()  { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed ()   { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers ()        { throw ExceptionNotImplemented(); }
};

/********************************************************************************/

class SystemInfoWindows : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores ();

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName ();

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal ()  { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed ()   { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers ()        { throw ExceptionNotImplemented(); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_SYSTEM_IMPL_SYSTEM_COMMON_HPP_ */
