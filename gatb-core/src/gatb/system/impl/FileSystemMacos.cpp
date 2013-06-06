/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#ifdef __APPLE__

#include <gatb/system/impl/FileSystemMacos.hpp>

#include <sys/resource.h>
#include <sys/statvfs.h>
#include <stdlib.h>
#include <dirent.h>

#include <string>
#include <sstream>

using namespace std;

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/

class MacosFile : public CommonFile
{
public:

    /** Constructor. */
	MacosFile (const char* path, const char* mode) : CommonFile(path,mode)  { }

    /** \copydoc IFile::tell */
    u_int64_t tell ()  { return (isOpen() ? ftello (_handle) : 0); }

	/** \copydoc IFile::seeko */
    int seeko (u_int64_t offset, int whence)  {  return (isOpen() ? fseeko (_handle, offset, whence) : -1);  }
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t FileSystemMacos::getMaxFilesNumber ()
{
    size_t result = 0;

    struct rlimit lim;

    if (getrlimit (RLIMIT_NOFILE, &lim) == 0)   {  result = lim.rlim_cur;  }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFile* FileSystemMacos::newFile (const Path& path, const char* mode)
{
    return new MacosFile (path.c_str(), mode);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFile* FileSystemMacos::newFile (const Path& dirpath, const Path& filename, const char* mode)
{
    /** We build the full file path. */
    stringstream ss;
    ss << dirpath << "/" << filename;

    /** We create the file handle. */
    return newFile (ss.str(), mode);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int FileSystemMacos::clearCache ()
{
    return ::system("purge");
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* __MACOS__ */
