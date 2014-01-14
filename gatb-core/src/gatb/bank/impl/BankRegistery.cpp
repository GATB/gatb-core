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

#include <gatb/bank/impl/BankRegistery.hpp>

#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/bank/impl/BankBinary.hpp>

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
void BankRegistery::registerFactory (const std::string& name, IBankFactory* instance)
{
    /** We look whether the factory is already registered. */
    std::map <std::string, IBankFactory*>::iterator it = _factories.find (name);
    if (it == _factories.end())
    {
        _factories[name] = instance;
        instance->use();
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
IBankFactory* BankRegistery::getFactory (const std::string& name)
{
    IBankFactory* result = 0;

    if (name.empty())  {  result = _factories[_defaultName]; }
    else               {  result = _factories[name];         }

    if (result == 0)  { throw system::Exception ("Bank factory not registered for name '%s'", name.c_str()); }

    /** We return the factory. */
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
BankRegistery::BankRegistery ()
{
    /** We register most known factories. */
    registerFactory ("fasta",  new BankFastaFactory());
    registerFactory ("binary", new BankBinaryFactory());

    /** We set the default one. */
    _defaultName = "fasta";
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankRegistery::~BankRegistery ()
{
    for (std::map <std::string, IBankFactory*>::iterator it = _factories.begin (); it != _factories.end(); it++)
    {
        (it->second)->forget ();
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
