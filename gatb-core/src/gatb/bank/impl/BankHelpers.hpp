/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
 *****************************************************************************/

/** \file BankHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Helpers for managing IBank objects
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Interface for reading genomic databases.
 */
class BankHelper
{
public:

    /** Singleton method.
     * \return the singleton instance. */
    static BankHelper& singleton()  { static BankHelper instance; return instance; }

    /** Convert one bank into another one.
     * \param[in] in  : the bank to be converted
     * \param[in] out : the converted bank
     * \param[in] progress : listener getting conversion progression information
     * */
    tools::misc::IProperties* convert (IBank& in, IBank& out, tools::dp::impl::IteratorListener* progress=0);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_ */
