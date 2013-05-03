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

#ifndef _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_
#define _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_

#include <gatb/system/api/types.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/********************************************************************************/

/** \brief Interface providing some general information about the system.
 */
class ISystemInfo
{
public:

    /** Returns the version of the library.
     * \return the version. */
    virtual std::string version () = 0;

    /** Returns the number of available cores.
     * \return the number of cores. */
    virtual size_t getNbCores () = 0;

    /** Returns the host name.
     * \return the host name. */
    virtual std::string getHostName () = 0;

    /** Get the size (in bytes) of the physical memory
     * \return the physical memory size */
    virtual u_int64_t getMemoryPhysicalTotal () = 0;

    /** Get the size (in bytes) of the used physical memory
     * \return the used physical memory size */
    virtual u_int64_t getMemoryPhysicalUsed () = 0;

    /** Get the size (in bytes) of the buffers memory
     * \return the buffers memory size */
    virtual u_int64_t getMemoryBuffers () = 0;

    /** Destructor. */
    virtual ~ISystemInfo ()  {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_ */
