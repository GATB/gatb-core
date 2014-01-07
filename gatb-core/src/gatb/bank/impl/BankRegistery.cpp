/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
