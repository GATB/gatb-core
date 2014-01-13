/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file ISystemInfo.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface providing information about the operating system
 */

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
    virtual std::string getVersion () const = 0;

    /** Returns the date of the library generation.
     * \return the generation date. */
    virtual std::string getBuildDate () const = 0;

    /** Returns the compilation options
     * \return the compilation options. */
    virtual std::string getBuildOptions () const = 0;

    /** Returns the operating system name
     * \return the os name. */
    virtual std::string getOsName () const = 0;

    /** Returns the number of available cores.
     * \return the number of cores. */
    virtual size_t getNbCores () const = 0;

    /** Returns the host name.
     * \return the host name. */
    virtual std::string getHostName () const = 0;

    /** Returns home directory.
     * \return the home directory uri. */
    virtual std::string getHomeDirectory () const = 0;

    /** Get the size (in bytes) of the physical memory
     * \return the physical memory size */
    virtual u_int64_t getMemoryPhysicalTotal () const = 0;

    /** Get the size (in bytes) of the used physical memory
     * \return the used physical memory size */
    virtual u_int64_t getMemoryPhysicalUsed () const = 0;

    /** Get the size (in bytes) of the buffers memory
     * \return the buffers memory size */
    virtual u_int64_t getMemoryBuffers () const = 0;

    /** Destructor. */
    virtual ~ISystemInfo ()  {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_ */
