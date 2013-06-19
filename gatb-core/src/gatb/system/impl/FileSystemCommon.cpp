/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/FileSystemCommon.hpp>

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
u_int64_t  FileSystemCommon::getAvailableSpace (const Path& path)
{
    struct statvfs buffer;

    statvfs (path.c_str(), &buffer);

    u_int64_t available = (buffer.f_bavail * buffer.f_bsize) / 1024;

    return available;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFileSystem::Path FileSystemCommon::getCurrentDirectory ()
{
    char path[1000];
    char* buffer = getcwd (path, sizeof(path));

    if (buffer == 0)  {  throw ExceptionErrno ("unable to get current directory");  }

    return path;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFileSystem::Path FileSystemCommon::getTemporaryDirectory ()
{
    const char* dir = 0;

         if ( (dir = getenv ("TMPDIR"))  != 0)  {  return dir;    }
    else if ( (dir = getenv ("TMP"))     != 0)  {  return dir;    }
    else if ( (dir = getenv ("TEMPDIR")) != 0)  {  return dir;    }
    else                                        {  return "/tmp"; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool FileSystemCommon::doesExist (const Path& path)
{
    FILE* fp = fopen (path.c_str(), "rb");

   if (fp != NULL)
   {
       fclose (fp);
       return true;
   }
   else
   {
       return false;
   }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
u_int64_t  FileSystemCommon::getSize (const Path& path)
{
    struct stat st;

    if (stat (path.c_str(), &st) == 0) return st.st_size;

    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void FileSystemCommon::iterate (const Path& path, void (*callback) (const Path&, void* data), void* data)
{
    DIR* dp = opendir (path.c_str());

    if (dp)
    {
        struct dirent* dirp = 0;

        while ( (dirp = readdir(dp)) != 0)
        {
            callback (dirp->d_name, data);
        }

        closedir (dp);
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
#if 0
IFile* FileSystemCommon::newFile (const Path& dirpath, const Path& filename, const char* mode)
{
    /** We build the full file path. */
    stringstream ss;
    ss << dirpath << "/" << filename;

    /** We create the file handle. */
    return newFile (ss.str(), mode);
}
#endif

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
