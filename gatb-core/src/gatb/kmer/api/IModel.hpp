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
#include <gatb/tools/misc/api/Data.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for kmer management. */
namespace kmer      {
/********************************************************************************/

/** \brief enumeration giving the way the kmers are computed. */
enum KmerMode
{
    /** Kmer from direct strand */
    KMER_DIRECT,
    /** Kmer as reverse complement */
    KMER_REVCOMP,
    /** Kmer as minimum between the direct and reverse complement */
    KMER_MINIMUM
};

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
     * \param[in] mode : mode telling how the kmer is computed
     * \return the kmer for the given nucleotides. */
    virtual kmer_type codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding) = 0;

    /** Compute the next right kmer given a current kmer and a nucleotide.
     * \param[in] val_seed : the current kmer as a starting point
     * \param[in] nucleotide : the next nucleotide
     * \param[in] mode : mode telling how the kmer is computed
     * \return the kmer on the right of the given kmer. */
    virtual kmer_type codeSeedRight (const kmer_type& val_seed, char nucleotide, tools::misc::Data::Encoding_e encoding) = 0;

    /** Build a vector of kmers from some source data (as a Data instance).
     * According to the implementation of the IModel interface, the kmers can be different (direct, revcomp, minimum)
     * \param[in] d : source data from which the kmers are built
     * \param[in] kmersBuffer : vector where to put the built kmers into.
     */
    virtual void build (tools::misc::Data& d, std::vector<kmer_type>& kmersBuffer) = 0;

    /** Destructor. */
    virtual ~IModel () {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMODEL_HPP_ */
