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

/** \file BankAlbum.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_ALBUM_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_ALBUM_HPP_

/********************************************************************************/

#include <gatb/bank/impl/BankComposite.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Interface for reading genomic databases.
 *
 * This class allows to define a list of banks URI in two ways:
 *  1) the album URI is a comma separated list of banks URI
 *  2) the album is a text file holding a list of banks URI
 *
 * Using a BankAlbum object is the same as using a list of banks defined by the banks.
 *
 * An album A is defined by N filepath Pi.
 * An album is valid if each Pi is present in filesystem
 *
 * If a filepath Pi is a simple filename (ie. not a fullpath), the basedir of the album URI
 * is prefixed to Pi. This feature avoids to write the fullpath in Pi, which is interesting
 * if the banks referred by Pi are moved; in such a case, one just has to move the album file
 * in the same location of the moved banks.
 *
 * BankAlbum is a composite bank.
 */
class BankAlbum : public BankComposite
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "album"; }

    /** Constructor.
     * \param[in] name : uri of the album. */
    BankAlbum (const std::string& name, bool deleteIfExists=false);

    /** \copydoc IBank::getId. */
    std::string getId ()  { return _name; }

    /** Add a bank to the album. */
    IBank* addBank (const std::string& bankUri);

    /** Add a bank to the album. */
    IBank* addBank (const std::string& directory, const std::string& bankName);

private:

    std::string _name;

    std::vector<std::string> _banksUri;

    system::IFile* getFile (const std::string& name, const char* mode=NULL);

    static bool isOnlyFilename (const std::string& path);

    /** */
    static bool isAlbumValid (const std::string& uri);

    friend class BankAlbumFactory;
};

/********************************************************************************/

/** \brief Factory for the BankAlbum class. */
class BankAlbumFactory : public IBankFactory
{
public:

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_ALBUM_HPP_ */
