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

/** \file Bank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of IBank interface with compressed format.
 */

#ifndef _GATB_CORE_BANK__IMPL_BANK_BINARY_HPP_
#define _GATB_CORE_BANK__IMPL_BANK_BINARY_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>

#include <vector>
#include <string>
#include <stdio.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IBank for compressed format
 */
class BankBinary : public IBank
{
public:

    /** Constructor.
     * \param[in] filename : uri of the bank. */
    BankBinary (const std::string& filename);

    /** Destructor. */
    ~BankBinary ();

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return new Iterator (*this); }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item);

    /** \copydoc IBank::flush */
    void flush ();

    /** \copydoc IBank::getSize */
    u_int64_t getSize ();

    /** \copydoc IBank::estimateNbSequences */
    u_int64_t estimateNbSequences ();

    /** \copydoc IBank::estimateNbSequences */
    size_t estimateMaxSequenceLength ()  { return 0; }

    /************************************************************/

    /** \brief Specific Iterator impl for BankBinary class
     */
    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:
        /** Constructor.
         * \param[in] ref : the associated iterable instance.
         */
        Iterator (BankBinary& ref);

        /** Destructor */
        ~Iterator ();

        /** \copydoc tools::dp::Iterator::first */
        void first();

        /** \copydoc tools::dp::Iterator::next */
        void next();

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _isDone; }

        /** \copydoc tools::dp::Iterator::item */
        Sequence& item ()     { return _sequence; }

        /** Estimation of the number of sequences. Used (by delegation) by the Bank class.
         * \return the sequences number estimation. */
        u_int64_t estimateNbSequences ();

    private:

        /** Reference to the underlying Iterable instance. */
        BankBinary&    _ref;

        /** We define a custom Sequence type in order to get dummy comments. */
        struct CustomSequence : bank::Sequence
        {
            /**  */
            const char* getComment () const
            {
                static char buffer[256];
                snprintf (buffer, sizeof(buffer), ">seq_%d_length_%ld", _idx, getDataSize());
                return buffer;
            }
            u_int32_t _idx;
        };

        /** Current item to be returned by the iterator. */
        CustomSequence _sequence;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        char* buffer;
        int   cpt_buffer;
        int   blocksize_toread;
        int   nseq_lues;

        unsigned int  read_write_buffer_size;

        FILE* binary_read_file;
    };

protected:

    /** URI of the bank. */
    std::string _filename;

    unsigned char* buffer;
    int            cpt_buffer;
    unsigned int   read_write_buffer_size;
    FILE*          binary_read_file;

    void open  (bool write);
    void close ();
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK__IMPL_BANK_BINARY_HPP_ */
