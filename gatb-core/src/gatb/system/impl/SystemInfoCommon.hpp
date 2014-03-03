/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include <gatb/system/api/IMemory.hpp>
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

    /** \copydoc ISystemInfo::getVersion */
    std::string getVersion () const  { return STR_LIBRARY_VERSION; }

    /** \copydoc ISystemInfo::getBuildDate */
    std::string getBuildDate () const { return STR_COMPILATION_DATE; }

    /** \copydoc ISystemInfo::getBuildCompiler */
    std::string getBuildCompiler () const  { return STR_COMPILER; }

    /** \copydoc ISystemInfo::buildOptions */
    std::string getBuildOptions () const { return STR_COMPILATION_FLAGS; }

    /** \copydoc ISystemInfo::getOsName */
    std::string getBuildSystem () const { return STR_OPERATING_SYSTEM; }

    /** \copydoc ISystemInfo::getHomeDirectory */
    std::string getHomeDirectory ()  const {  return getenv("HOME") ? getenv("HOME") : ".";  }

    /** \copydoc ISystemInfo::getMemoryProject */
    u_int64_t getMemoryPhysicalFree () const  { return getMemoryPhysicalTotal()-getMemoryPhysicalUsed(); }

    /** \copydoc ISystemInfo::getMemoryProject */
    u_int64_t getMemoryProject () const  {  return std::min (getMemoryPhysicalFree() / (2*MBYTE), (u_int64_t)(5*1024)); }
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
    u_int64_t getMemoryPhysicalTotal () const;

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed () const;

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
