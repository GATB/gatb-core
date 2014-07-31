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
#include <list>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Class providing factories for building IBank objects.
 *
 * This class registers IBankFactory instances, which allows to create IBank objects
 * through a call to 'createBank'.
 *
 * Today, the following factories are registered:
 *  1) BankAlbumFactory
 *  2) BankFastaFactory
 *  3) BankBinaryFactory
 *
 * During a call to 'createBank', each factory is tried (in the order of registration)
 * until a correct IBank object is returned; if no valid IBank is found, an exception
 * is thrown.
 *
 * This class is a Singleton (private constructor) and must be used only this way.
 */
class BankRegistery
{
public:

    /** Singleton instance. */
    static BankRegistery& singleton()  { static BankRegistery instance; return instance; }

    /** Register a new factory, associated with a name.
     * \param[in] name : name of the factory
     * \param[in] instance : IBank factory */
    void registerFactory (const std::string& name, IBankFactory* instance);

    /** Get a bank handle.
     * \param[in] uri : uri of the bank.
     * \return the bank handle. */
    IBank* createBank (const std::string& uri);

    /** Get the type of the bank as a string
     * \param[in] uri : uri of the bank.
     * \return the bank type as a string. */
    std::string getType (const std::string& uri);

    /** Get a factory for a given name. */
    IBankFactory* getFactory (const std::string& name);

private:

    /** Private due to singleton method. */
    BankRegistery ();

    /** Destructor. */
    ~BankRegistery ();

    /** The order of registration is important, so we can't rely on a map and have to use
     * a list of Entry (equivalent to <key,value> of a map). */
    struct Entry
    {
        Entry (const std::string& name, IBankFactory* factory) : name(name), factory(factory) {}
        std::string   name;
        IBankFactory* factory;
    };

    std::list<Entry> _factories;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_ */
