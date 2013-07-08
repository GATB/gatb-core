/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file IBank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Definition of what a genomic sequence is.
 */

#ifndef _GATB_CORE_BANK_SEQUENCE_HPP_
#define _GATB_CORE_BANK_SEQUENCE_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/Data.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
/********************************************************************************/

/** \brief Interface for what we need to read genomic databases.
 */
struct Sequence
{
    /** Constructor. */
    Sequence (tools::misc::Data::Encoding_e encoding = tools::misc::Data::ASCII) : _data(encoding)  {}

    /** Destructor. */
    virtual ~Sequence ()  { }

    /** \return description of the sequence */
    virtual const std::string& getComment ()  const  { return _comment; }

    /** \return the data as a Data structure. */
    virtual tools::misc::Data& getData () { return _data; }

    /** \return the data as a Data structure. */
    virtual u_int64_t getSeqNum () { return _seqNum; }
    
    /** \return buffer holding the sequence residues. */
    virtual char* getDataBuffer ()  const { return _data.getBuffer(); }

    /** \return number of residues of the sequence. */
    virtual size_t getDataSize () const  { return _data.size(); }

    /** \return format of the data. */
    virtual tools::misc::Data::Encoding_e getDataEncoding () const  { return _data.getEncoding(); }

    void setDataRef (tools::misc::Data* ref, int offset, int length)
    {
        _data.setRef (ref, offset, length);
    }

    std::string _comment;

    u_int64_t _seqNum;

private:
    tools::misc::Data _data;
    
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_SEQUENCE_HPP_ */
