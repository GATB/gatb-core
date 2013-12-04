/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
    ~BankConverterAlgorithm ();

    /** */
    void execute ();

    /** */
    IBank* getResult ()  { return _bankOutput; }

private:

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

