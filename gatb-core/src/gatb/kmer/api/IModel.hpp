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

/** \file IModel.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for the kmer model
 */

#ifndef _GATB_CORE_KMER_IMODEL_HPP_
#define _GATB_CORE_KMER_IMODEL_HPP_

/********************************************************************************/

#include <gatb/bank/api/IAlphabet.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for kmer management. */
namespace kmer      {
/********************************************************************************/

/** \brief High level interface for dealing with kmer
 *
 * The IModel provides a mean to do several actions on kmer such as:
 *
 * 1) get the size of the kmer
 * 2) get the underlying alphabet
 * 3) getting the next kmer from one kmer and the next letter (on the right or on the left)
 */
template <typename kmer_type> class IModel
{
public:

    /** \return the alphabet. */
    virtual bank::IAlphabet&  getAlphabet() = 0;

    /** \return the span (ie the kmer size) for this model. */
    virtual size_t getSpan () = 0;

    /** \return the memory size of a kmer. */
    virtual size_t getMemorySize () = 0;

    /** Compute the kmer given some nucleotide data.
     *  Note that we don't check if we have enough nucleotides in the provided data.
     * \param[in] seq : the sequence
     * \param[in] revcomp : reverse code if true, direct otherwise
     * \return the kmer for the given nucleotides. */
    virtual kmer_type codeSeed (const char* seq, bool revcomp) = 0;

    /** Compute the next right kmer given a current kmer and a nucleotide.
     * \param[in] val_seed : the current kmer as a starting point
     * \param[in] nucleotide : the next nucleotide
     * \param[in] revcomp : reverse code if true, direct otherwise
     * \return the kmer on the right of the given kmer. */
    virtual kmer_type codeSeedRight (const kmer_type& val_seed, char nucleotide, bool revcomp) = 0;

    /** Destructor. */
    virtual ~IModel () {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMODEL_HPP_ */
