/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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

#include <gatb/bank/impl/BankAlbum.hpp>
#include <gatb/bank/impl/BankRegistery.hpp>
#include <gatb/system/impl/System.hpp>

#include <fstream>

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankAlbum::BankAlbum (const std::string& name) : _name(name)
{
    /** We get a handle on the file given its name. */
    system::IFile* file = getFile (name);

    /** We check that the provided name exists in filesystem. */
    if (file != 0)
    {
        char buffer[256];

        while (file->isEOF() == false)
        {
            /** We init the buffer. */
            *buffer = 0;

            /** We get the current line. */
            int len = file->gets (buffer, sizeof(buffer));

            if (len > 1)
            {
                /** We remove the end of line. */
                if (buffer[len-1] == '\n')  {  buffer[len-1] = 0; }

                string bankUri = buffer;

                /** We check whether it is a mere file name or there is also a directory name. */
                if (isOnlyFilename(buffer) == true)
                {
                    /** We add the path name of the album file. */
                    bankUri = System::file().getDirectory (name) + "/" + bankUri;
                }

                /** We add a new bank. */
                addBank (BankRegistery::singleton().getFactory()->createBank(bankUri));
            }
        }

        delete file;
    }
    else
    {
        throw Exception ("Unable to use file '%s'", name.c_str());
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
BankAlbum::~BankAlbum ()
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankAlbum::add (const std::string& bankUri)
{
    DEBUG (("BankAlbum::add '%s'\n", bankUri.c_str() ));

    /** We add the uri into the album file. */
    system::IFile* file = getFile (_name);

    if (file != 0)
    {
        /** We write the uri in the file. */
        file->print ("%s\n", bankUri.c_str());

        /** We add a new bank. */
        addBank (BankRegistery::singleton().getFactory()->createBank(bankUri));

        /** Cleanup. */
        delete file;
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
system::IFile* BankAlbum::getFile (const std::string& name)
{
    /** We check whether the file already exists or not. */
    const char* mode = System::file().doesExist(name) == true ? "r+" : "w+";

    /** We create a handle on the file. */
    return System::file().newFile (name, mode);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankAlbum::isOnlyFilename (const std::string& path)
{
    if (path.empty())  { throw Exception ("Bad '%s' path in isOnlyFilename", path.c_str());  }

    /** It may not be bullet proof... */
    return (path[0]!='/' && path[0]!='.');
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
