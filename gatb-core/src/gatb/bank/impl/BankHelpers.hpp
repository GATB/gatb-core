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

/** \file BankHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Helpers for managing IBank objects
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>
#include <gatb/bank/impl/BankRegistery.hpp>
#include <gatb/bank/impl/AbstractBank.hpp>
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
    tools::misc::IProperties* convert (IBank& in, IBank& out, tools::dp::IteratorListener* progress=0);
};

/********************************************************************************/

class BankDelegate : public AbstractBank
{
public:

    /** Constructor.
     * \param[in] ref : referred bank.
     */
    BankDelegate (IBank* ref) : _ref(0)  { setRef(ref); }

    /** Destructor. */
    ~BankDelegate () { setRef(0); }

    /** */
    std::string getId ()  { return _ref->getId(); }

    /** \copydoc tools::collections::Iterable::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return _ref->iterator(); }

    /**  */
    int64_t getNbItems () { return _ref->getNbItems(); }

    /** \copydoc tools::collections::Bag */
    void insert (const Sequence& item)   { _ref->insert (item); }

    /** */
    void flush ()  {  _ref->flush ();  }

    /** Return the size of the bank (comments + data)
     *
     * The returned value may be an approximation in some case. For instance, if we use
     * a zipped bank, an implementation may be not able to give accurate answer to the
     * size of the original file.
     *
     * \return the bank size.*/
    u_int64_t getSize ()  { return _ref->getSize(); }

    /** Give an estimation of sequences information in the bank:
     *      - sequences number
     *      - sequences size (in bytes)
     *      - max size size (in bytes)
     * \return the sequences number estimation. */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)  {  _ref->estimate (number, totalSize, maxSize);  }

    /** Shortcut to 'estimate' method.
     * \return estimation of the number of sequences */
    int64_t estimateNbItems () { return _ref->estimateNbItems(); }

    /** Shortcut to 'estimate' method.
     * \return estimation of the size of sequences */
    u_int64_t estimateSequencesSize ()  { return _ref->estimateSequencesSize(); }

    /** \return the number of sequences read from the bank for computing estimated information */
    u_int64_t getEstimateThreshold ()  { return _ref->getEstimateThreshold(); }

    /** Set the number of sequences read from the bank for computing estimated information
     * \param[in] nbSeq : the number of sequences to be read.*/
    void setEstimateThreshold (u_int64_t nbSeq)  { _ref->setEstimateThreshold(nbSeq); }

protected:

    IBank* _ref;
    void setRef (IBank* ref)  { SP_SETATTR(ref); }
};

/********************************************************************************/

template<typename Filter> class BankFiltered : public BankDelegate
{
public:

    /** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
    BankFiltered (IBank* ref, const Filter& filter) : BankDelegate (ref), _filter(filter)  {}

    /** \copydoc tools::collections::Iterable::iterator */
    tools::dp::Iterator<Sequence>* iterator ()
    {
        return new tools::dp::impl::FilterIterator<Sequence,Filter> (_ref->iterator (), _filter);
    }

private:

    Filter _filter;
};

/********************************************************************************/

template<typename Filter> class BankFilteredFactory : public IBankFactory
{
public:

    BankFilteredFactory (const std::string& delegateFormat, const Filter& filter) : _format(delegateFormat), _filter(filter)  {}

    IBank* createBank (const std::string& uri)
    {
        /** We create the reference bank. */
        IBank* ref = BankRegistery::singleton().getFactory(_format)->createBank (uri);

        /** We encapsulate with a filtered bank. */
        return new BankFiltered<Filter> (ref, _filter);
    }

private:

    std::string _format;
    Filter      _filter;
};


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_ */
