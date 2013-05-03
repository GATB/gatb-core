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

/** \file Alphabet.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation for genomic alphabets
 */

#ifndef _GATB_CORE_BANK__IMPL_ALPHABET_HPP_
#define _GATB_CORE_BANK__IMPL_ALPHABET_HPP_

/********************************************************************************/

#include <gatb/bank/api/IAlphabet.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** Define the kind of the underlying alphabet. For instance, it could be the kind of sequences
 * (protein, ADN...) read from a FASTA file
 */
class AlphabetNucleic : public IAlphabet
{
public:

    /** Singleton method.
     * \return the singleton instance. */
    static IAlphabet& singleton()  { static AlphabetNucleic instance; return instance; }

    /** \copydoc IAlphabet::getKind */
    Kind_e getKind ()  { return NUCLEIC_ACID; }

    /** \copydoc IAlphabet::getLetters */
    std::string  getLetters ()  { return std::string ("ACTG"); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK__IMPL_ALPHABET_HPP_ */
