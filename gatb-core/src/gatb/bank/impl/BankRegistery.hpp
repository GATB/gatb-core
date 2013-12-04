/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
