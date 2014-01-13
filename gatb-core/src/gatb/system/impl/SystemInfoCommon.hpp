/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
#include <gatb/system/api/config.hpp>

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
    std::string getVersion () const  { return STR_LIBRARY_VERSION; }

    /** \copydoc ISystemInfo::buildDate */
    std::string getBuildDate () const { return STR_COMPILATION_DATE; }

    /** \copydoc ISystemInfo::buildOptions */
    std::string getBuildOptions () const { return STR_COMPILATION_FLAGS; }

    /** \copydoc ISystemInfo::getOsName */
    std::string getOsName () const { return STR_OPERATING_SYSTEM; }

    /** \copydoc ISystemInfo::getHomeDirectory */
    std::string getHomeDirectory ()  const {  return getenv("HOME") ? getenv("HOME") : ".";  }
};

/********************************************************************************/

class SystemInfoLinux : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores () const ;

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName () const ;

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal () const ;

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed () const ;

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers () const ;
};

/********************************************************************************/

class SystemInfoMacos : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores () const ;

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName () const ;

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal () const  { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed () const   { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers () const        { throw ExceptionNotImplemented(); }
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
    u_int64_t getMemoryPhysicalTotal () const  { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed ()  const   { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers ()       const  { throw ExceptionNotImplemented(); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_SYSTEM_IMPL_SYSTEM_COMMON_HPP_ */
