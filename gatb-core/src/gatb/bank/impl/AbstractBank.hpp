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

/** \file AbstractBank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Abstract implementation of the IBank interface
 */

#ifndef _GATB_CORE_BANK_IMPL_ABSTRACT_BANK_HPP_
#define _GATB_CORE_BANK_IMPL_ABSTRACT_BANK_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Abstract implementation of IBank for factorizing common behavior.
 */
class AbstractBank : public IBank
{
public:

    /** Constructor. */
    AbstractBank () : _estimateThreshold(5000) {}

    /** \copydoc IBank::estimateNbSequences */
    u_int64_t estimateNbSequences ()
    {
        u_int64_t number, totalSize, maxSize;    estimate (number, totalSize, maxSize);  return number;
    }

    /** \copydoc IBank::estimateSequencesSize */
    u_int64_t estimateSequencesSize ()
    {
        u_int64_t number, totalSize, maxSize;    estimate (number, totalSize, maxSize);  return totalSize;
    }

    /** \copydoc IBank::getEstimateThreshold */
    u_int64_t getEstimateThreshold ()  { return _estimateThreshold; }

    /** \copydoc IBank::setEstimateThreshold */
    void setEstimateThreshold (u_int64_t nbSeq) { _estimateThreshold = nbSeq; }

private:

    u_int64_t _estimateThreshold;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_ABSTRACT_BANK_HPP_ */
