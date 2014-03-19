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

/** \file BankFasta.hpp
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
 */
class BankAlbum : public BankComposite
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "album"; }

    /** Constructor.
     * \param[in] name : uri of the album. */
    BankAlbum (const std::string& name, bool deleteIfExists=false);

    /** Destructor */
    virtual ~BankAlbum ();

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

    bool isOnlyFilename (const std::string& path);
};

/********************************************************************************/

/** */
class BankAlbumFactory : public IBankFactory
{
public:

    IBank* createBank (const std::string& uri) { return new BankAlbum (uri); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_ALBUM_HPP_ */
