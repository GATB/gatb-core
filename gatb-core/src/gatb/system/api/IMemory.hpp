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

/** \file IMemory.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of an abstraction for dynamic memory allocation.
 *
 * We define abstraction for malloc,realloc... functions. Note that failures
 * during execution of the methods should throw a specific exception.
 */

#ifndef _GATB_CORE_SYSTEM_IMEMORY_HPP_
#define _GATB_CORE_SYSTEM_IMEMORY_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/********************************************************************************/

/** \brief enumeration defining memory units */
enum MemoryUnit
{
    KBYTE = (1ULL << 10),
    MBYTE = (1ULL << 20),
    GBYTE = (1ULL << 30)
};

/********************************************************************************/

/** \brief Interface providing methods for dynamic allocation.
 *
 *  This interface provides most common methods for creating/deleting dynamic allocation buffers.
 */
class IMemoryAllocator
{
public:

    /************************************************************/
    /** Alias for the size of memory blocks handled by the interface.  */
    typedef u_int32_t BlockSize_t;

    /** Alias for the size of memory handled by the interface.  */
    typedef u_int64_t TotalSize_t;

    /************************************************************/
    /** See malloc documentation. */
    virtual void* malloc  (BlockSize_t size) = 0;

    /** See calloc documentation. */
    virtual void* calloc  (size_t nmemb, BlockSize_t size) = 0;

    /** See realloc documentation. */
    virtual void* realloc (void *ptr, BlockSize_t size) = 0;

    /** See free documentation. */
    virtual void  free    (void *ptr) = 0;

    /************************************************************/
    /** Destructor. */
    virtual ~IMemoryAllocator () {}
};

/********************************************************************************/

/** \brief Interface providing methods for manipulating memory blocks
 *
 *  This interface provides most common methods for setting/copying/comparing buffers.
 */
class IMemoryOperations
{
public:

    /** Same as memset from <string.h> */
    virtual void* memset (void* s, int c, size_t n) = 0;

    /** Same as memcpy from <string.h> */
    virtual void* memcpy (void* dest, const void* src, size_t n) = 0;

    /** Same as memcpy from <string.h> */
    virtual int memcmp (const void* s1, const void* s2, size_t n) = 0;

    /************************************************************/
    /** Destructor. */
    virtual ~IMemoryOperations () {}
};

/********************************************************************************/

/** \brief Interface providing methods for dynamic allocation and memory statistics
 *
 *  This interface provides most common methods for creating/deleting dynamic allocation buffers.
 *
 *  It also provides information about its inner state (current memory usage, maximum reached).
 *
 *  Specific implementations could be:
 *      - system allocator (ie. functions from stdlib.h)
 *      - allocator with a maximum allowed total size
 *      - etc...
 */
class IMemory : public IMemoryAllocator, public IMemoryOperations
{
public:

    /** Get the number of currently allocated memory blocks.
     * \return the number of allocated blocks. */
    virtual size_t getNbBlocks () = 0;

    /** Get the memory usage by the current process.
     * \return the memory usage (in bytes). */
    virtual TotalSize_t getCurrentUsage () = 0;

    /** Get the maximum reached size.
     * \return the maximum used memory (in bytes). */
    virtual TotalSize_t getMaximumUsage () = 0;

    /************************************************************/
    /** Destructor. */
    virtual ~IMemory () {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IMEMORY_HPP_ */
