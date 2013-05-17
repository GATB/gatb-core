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

#include <gatb/bank/impl/AbstractBank.hpp>

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
 *
 * - a binary file is a list of blocks
 *    - a block is:
 *       - one block size (on 4 bytes)
 *       - a list of sequences
 *          - a sequence is:
 *             - a sequence length (on 4 bytes)
 *             - the nucleotides of the sequences (4 nucleotides encoded in 1 byte)
 * - number of sequences (on 4 bytes)
 */
class BankBinary : public AbstractBank
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

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

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
        virtual ~Iterator ();

        /** \copydoc tools::dp::Iterator::first */
        void first();

        /** \copydoc tools::dp::Iterator::next */
        void next();

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _isDone; }

        /** \copydoc tools::dp::Iterator::item */
        Sequence& item ()     { return *_item; }

        /** Estimation of the sequences information. */
        void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    private:

        /** Reference to the underlying Iterable instance. */
        BankBinary&    _ref;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        /** Block buffer read from file. */
        tools::misc::Data* _bufferData;
        void setBufferData (tools::misc::Data* bufferData)  { SP_SETATTR(bufferData); }

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
