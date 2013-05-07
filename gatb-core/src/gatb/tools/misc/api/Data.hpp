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

/** \file Range.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Data structure
 */

#ifndef _GATB_CORE_TOOLS_MISC_DATA_HPP_
#define _GATB_CORE_TOOLS_MISC_DATA_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief Definition of a data chunk
 *
 * A data is defined by:
 *      - an encoding format
 *      - a buffer holding the actual data
 *      - the size of the data
 */
struct Data
{
    /** Define how data is encoded. */
    enum Encoding_e
    {
        /** data encoded as ASCII codes (so one byte per data unit) */
        ASCII,
        /** one byte per data as integer value (for instance: A=0, C=1, T=2, G=3) */
        INTEGER,
        /** 4 nucleotides compressed in one byte */
        BINARY
    };

    /** Constructor.
     * \param[in] buf : the buffer holding the actual data
     * \param[in] len : the size of the data
     * \param[in] encode : the encoding scheme of the buffer
     */
    Data (char* buf, size_t len, Encoding_e encode)  : buffer(buf), size(len), encoding(encode)  {}

    /** Default constructor. */
    Data ()  : buffer(0), size(0), encoding(ASCII)  {}

    /** \return buffer holding the actual data. */
    char* getBuffer () const  { return buffer; }

    /** \return buffer size (in bytes). */
    size_t getSize ()  const  { return size; }

    /** \return format of the data. */
    Encoding_e getEncoding ()  const  { return encoding; }

    /** Operator overload for accessing a specific data in the buffer.
     * \param[in] idx : index of the data to be retrieved
     * \return the ith data
     */
    char& operator[] (size_t idx)  { return buffer[idx]; }

    char*       buffer;
    int         size;
    Encoding_e  encoding;
};

/********************************************************************************/

struct DataConverter
{
    static void convert (Data& in, Data& out)
    {
        size_t nchar = (in.getSize()+3)/4;
        size_t j=0;
        for (size_t i=0; i<nchar; i++)
        {
            char fournt = in[i];
            out[j+3] = fournt & 3; fournt = fournt >> 2;
            out[j+2] = fournt & 3; fournt = fournt >> 2;
            out[j+1] = fournt & 3; fournt = fournt >> 2;
            out[j+0] = fournt & 3;
            j+=4;
        }

        out.encoding = Data::ASCII;
        out.size     = in.size;
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_DATA_HPP_ */
