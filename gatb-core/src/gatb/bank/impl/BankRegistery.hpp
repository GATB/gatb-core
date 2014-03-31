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

/** \file BankRegistery.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>

#include <string>
#include <map>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief TBD
 */
class BankRegistery
{
public:

    /** Singleton instance. */
    static BankRegistery& singleton()  { static BankRegistery instance; return instance; }

    /** */
    void registerFactory (const std::string& name, IBankFactory* instance);

    /** */
    IBankFactory* getFactory (const std::string& name = "");

    /** Get a bank handle.
     * \param[in] uri : uri of the bank.
     * \return the bank handle. */
    IBank* createBank (const std::string& uri)  {  return getFactory()->createBank(uri);  }

private:

    /** Private due to singleton method. */
    BankRegistery ();

    /** */
    ~BankRegistery ();

    std::map <std::string, IBankFactory*> _factories;

    std::string _defaultName;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_ */
