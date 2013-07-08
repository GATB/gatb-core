/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file IBank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK_IBANK_HPP_
#define _GATB_CORE_BANK_IBANK_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/bank/api/Sequence.hpp>

/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writting
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/********************************************************************************/

/** \brief Interface for what we need to read genomic databases.
 *
 * The IBank interface is a essentially a short name for an Iterable over Sequence
 * instances.
 */
class IBank : public tools::collections::Iterable<Sequence>, public tools::collections::Bag<Sequence>
{
public:

    /** \copydoc tools::collections::Iterable::iterator */
    virtual tools::dp::Iterator<Sequence>* iterator () = 0;

    /** \copydoc tools::collections::Bag */
    virtual void insert (const Sequence& item) = 0;

    /** Return the size of the bank (comments + data)
     *
     * The returned value may be an approximation in some case. For instance, if we use
     * a zipped bank, an implementation may be not able to give accurate answer to the
     * size of the original file.
     *
     * \return the bank size.*/
    virtual u_int64_t getSize () = 0;

    /** Give an estimation of sequences information in the bank:
     *      - sequences number
     *      - sequences size (in bytes)
     *      - max size size (in bytes)
     * \return the sequences number estimation. */
    virtual void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize) = 0;

    /** Shortcut to 'estimate' method.
     * \return estimation of the number of sequences */
    virtual u_int64_t estimateNbSequences () = 0;

    /** Shortcut to 'estimate' method.
     * \return estimation of the size of sequences */
    virtual u_int64_t estimateSequencesSize () = 0;

    /** \return the number of sequences read from the bank for computing estimated information */
    virtual u_int64_t getEstimateThreshold () = 0;

    /** Set the number of sequences read from the bank for computing estimated information
     * \param[in] nbSeq : the number of sequences to be read.*/
    virtual void setEstimateThreshold (u_int64_t nbSeq) = 0;
};

/********************************************************************************/

/** \brief Factory for IBank.
 */

class IBankFactory : public tools::dp::SmartPointer
{
public:

    virtual IBank* createBank (const std::string& uri) = 0;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IBANK_HPP_ */
