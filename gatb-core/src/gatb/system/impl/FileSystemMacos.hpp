/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file FileSystemMacos.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation for MacOs.
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_FILE_SYSTEM_MACOS_HPP_
#define _GATB_CORE_SYSTEM_IMPL_FILE_SYSTEM_MACOS_HPP_

/********************************************************************************/

#include <gatb/system/impl/FileSystemCommon.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {
/********************************************************************************/

/** \brief default implementation
 */
class FileSystemMacos : public FileSystemCommon
{
public:

	/** \copydoc IFileSystem::getMaxFilesNumber */
    size_t getMaxFilesNumber ();

    /** \copydoc IFileSystem::newFile */
    IFile* newFile (const Path& path, const char* mode);

    /** \copydoc IFileSystem::newFile(const Path&, const Path&, const char*) */
    IFile* newFile (const Path& dirpath, const Path& filename, const char* mode);

    /** \copydoc IFileSystem::clearCache */
    int clearCache ();
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IMPL_FILE_SYSTEM_MACOS_HPP_ */
