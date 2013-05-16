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

/** \file ModelAbstract.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Basic implementation of the IModel interface
 */

#ifndef _GATB_CORE_KMER_IMPL_MODEL_ABSTRACT_HPP_
#define _GATB_CORE_KMER_IMPL_MODEL_ABSTRACT_HPP_

/********************************************************************************/

#include <gatb/kmer/api/IModel.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/misc/api/Data.hpp>

/** Note: we need to include here all potential integer type used by the template
 * parameter of the ModelAbstract class. */
#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>

/** External static data. */
extern const unsigned char comp_NT    [];

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** \brief Partial implementation of the IModel interface.
 *
 * This implementation factorizes some code for its children classes.
 *
 * In particular, it provides a generic Iterator implementation relying on the IModel::build
 * method called by the setData method.
 *
 * It also provides a generic 'codeSeed' method that computes a kmer from some buffer for the
 * given mode (direct kmer, revcomp or the minimum of both). This method is likely to be used
 * by children classes for implementing the methods not implemented here.
 */
template <typename kmer_type> class ModelAbstract : public IModel<kmer_type>
{
public:

    /** Constructor.
     * \param[in] span : size of the kmers for this model
     * \param[in] alphabet : underlying alphabet
     */
    ModelAbstract (size_t span, bank::IAlphabet& alphabet) : _sizeKmer(span), _alphabet(alphabet)
    {
        /** We compute the mask of the kmer. Useful for computing kmers in a recursive way. */
        _kmerMask = (((kmer_type)1) << (_sizeKmer*2))-1;

        size_t shift = 2*(_sizeKmer-1);

        /** The _revcompTable is a shortcut used while computing revcomp recursively. */

        /** Important: don't forget the kmer_type cast, otherwise the result in only on 32 bits. */
        for (size_t i=0; i<4; i++)  {  _revcompTable[i] = ((kmer_type)comp_NT[i]) << shift;  }
    }

    /** \copydoc IModel::getAlphabet */
    bank::IAlphabet&  getAlphabet() { return _alphabet; }

    /** \copydoc IModel::getSpan */
    size_t getSpan () { return _sizeKmer; }

    /** \copydoc IModel::getMemorySize */
    size_t getMemorySize ()  { return sizeof (kmer_type); }

    /************************************************************/
    /** \brief Specific Iterator impl for Model class */
    class Iterator : public tools::dp::impl::VectorIterator<kmer_type>
    {
    public:
        /** Constructor.
         * \param[in] ref : the associated model instance.
         */
        Iterator (ModelAbstract& ref)  : _ref(ref)   {}

        /** Set the data to be iterated.
         * \param[in] data : the data as information source for the iterator
         */
        void setData (tools::misc::Data& d)
        {
            /** We fill the vector with the items to be iterated. */
            _ref.build (d, this->_items);

            /** We set the vector size. */
            this->_nb = this->_items.size();
        }

    private:
        /** Reference on the underlying model; called for its 'build' method. */
        ModelAbstract& _ref;
    };

protected:

    /** Size of a kmer for this model. */
    size_t           _sizeKmer;

    /** Alphabet telling what kind of letters are supposed to be used.
     * Not really used so far, but put here for the future if we will also deal with
     * amino acids alphabets for implementing other algorithms than assembly. */
    bank::IAlphabet& _alphabet;

    /** Mask for the kmer. Used for computing recursively kmers. */
    kmer_type  _kmerMask;

    /** \return the kmer mask for the model, ie for the model kmer size. */
    kmer_type getMask () { return _kmerMask; }

    /** Shortcut for easing/speeding up the recursive revcomp computation. */
    kmer_type _revcompTable[4];

    /** Functor that returns the provided nucleotide.
     * \return the same value as input. */
    static int NTidentity(char nt)  { return nt; }

    /** Transform a nucleotide in ASCII form into an integer form as:
     *     - A=0
     *     - C=1
     *     - T=2
     *     - G=3
     * \return the translated nucleotide */
    static int NT2int(char nt)  {  return (nt>>1)&3;  }

    /** Generic function that computes a kmer from a sequence.
     * \param[in] seq : source data from which the kmer is computed.
     * \param[in] encoding : encoding of the source data
     * \param[in] mode : tells how to compute the kmer.
     * \return the computed kmer.
     */
    kmer_type codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding, KmerMode mode)
    {
        kmer_type direct = 0;

        if (encoding == tools::misc::Data::ASCII)  { for (size_t i=0; i<this->_sizeKmer; ++i)  {  direct = direct*4 + NT2int(seq[i]);  }  }
        else                                       { for (size_t i=0; i<this->_sizeKmer; ++i)  {  direct = direct*4 +       (seq[i]);  }  }

        if (mode == KMER_DIRECT)  { return direct; }

        kmer_type rev = core::tools::math::revcomp (direct, _sizeKmer);

        if (mode == KMER_REVCOMP)  { return rev; }

        return std::min (direct, rev);
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MODEL_ABSTRACT_HPP_ */
