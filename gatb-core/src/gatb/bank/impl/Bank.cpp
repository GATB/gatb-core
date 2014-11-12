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

#include <gatb/bank/impl/Bank.hpp>

#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankAlbum.hpp>

#include <gatb/system/api/Exception.hpp>

using namespace std;

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
Bank::Bank ()
{
    /** We register most known factories. */
    registerFactory ("album",  new BankAlbumFactory());
    registerFactory ("fasta",  new BankFastaFactory());
    registerFactory ("binary", new BankBinaryFactory());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bank::~Bank ()
{
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        (it->factory)->forget ();
    }
    _factories.clear();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bank::registerFactory (const std::string& name, IBankFactory* instance)
{
    /** We look whether the factory is already registered. */
    IBankFactory* factory = getFactory (name);
    if (factory == 0)
    {
        _factories.push_back (Entry (name, instance));
        instance->use();
    }
    else
    {
        throw string("Bank factory ") + name + string (" already registered");
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
IBankFactory* Bank::getFactory (const std::string& name)
{
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        if (it->name == name)  { return it->factory; }
    }
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
IBank* Bank::createBank (const std::string& uri)
{
    DEBUG (("Bank::createBank : %s \n", uri.c_str()));

    IBank* result = 0;
    for (list<Entry>::iterator it = _factories.begin(); result==0 && it != _factories.end(); it++)
    {
        result = it->factory->createBank(uri);
        DEBUG (("   factory '%s' => result=%p \n", it->name.c_str(), result ));
    }

    if (result == 0) { throw string("Bad bank creation from bank registery"); }

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
std::string Bank::getType (const std::string& uri)
{
    string result = "unknown";

    /** We try to create the bank; if a bank is valid, then we have the factory name. */
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        IBank* bank = it->factory->createBank(uri);
        if (bank != 0)
        {
            result = it->name;
            delete bank;
            break;
        }
    }

    return result;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
