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

/** \file BankConverterAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bank conversion from one IBank to another IBank
 */

#ifndef _BANK_CONVERTER_ALGORITHM_HPP_
#define _BANK_CONVERTER_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

class BankConverterAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** */
    BankConverterAlgorithm (IBank* bank, size_t kmerSize, const std::string& outputUri);

    /** */
    BankConverterAlgorithm (tools::storage::impl::Storage& storage);

    /** */
    ~BankConverterAlgorithm ();

    /** */
    void execute ();

    /** */
    IBank* getResult ()  { return _bankOutput; }

private:

    tools::misc::BankConvertKind _kind;

    IBank*      _bankInput;
    void setBankInput (IBank* bankInput)  { SP_SETATTR(bankInput); }

    IBank*      _bankOutput;
    void setBankOutput (IBank* bankOutput)  { SP_SETATTR(bankOutput); }

    std::string _outputUri;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _BANK_CONVERTER_ALGORITHM_HPP_ */

