/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Model.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Basic implementation of the IModel interface
 */

#ifndef _GATB_CORE_KMER_IMPL_MODEL_HPP_
#define _GATB_CORE_KMER_IMPL_MODEL_HPP_

/********************************************************************************/

#include <gatb/system/api/Exception.hpp>
#include <gatb/kmer/api/IModel.hpp>
#include <gatb/kmer/impl/ModelAbstract.hpp>

#include <gatb/bank/api/IAlphabet.hpp>
#include <gatb/bank/impl/Alphabet.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/Data.hpp>

#include <vector>
#include <algorithm>

#include <iostream>

extern const char bin2NT[] ;
extern const char binrev[] ;

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

extern const unsigned char revcomp_4NT[];
extern const unsigned char comp_NT    [];

    
/** \brief TO BE DONE
 */
template <typename kmer_type> class Model : public ModelAbstract<kmer_type>
{
public:

    /** Constructor.
     * \param[in] span : size of the kmers for this model
     * \param[in] alphabet : underlying alphabet
     */
    Model (size_t span, bank::IAlphabet& alphabet= bank::impl::AlphabetNucleic::singleton()) : ModelAbstract<kmer_type> (span, alphabet)  {}

    /** \copydoc IModel::codeSeed */
    kmer_type codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding)
    {
        return ModelAbstract<kmer_type>::codeSeed (seq, encoding, KMER_MINIMUM);
    }

    /** \copydoc IModel::codeSeedRight */
    kmer_type codeSeedRight (const kmer_type& val_seed, char nucleotide, tools::misc::Data::Encoding_e encoding,KmerMode mode = KMER_MINIMUM)
    {
        return ModelAbstract<kmer_type>::codeSeedRight (val_seed, nucleotide, encoding, mode);
    }

    /** */
    kmer_type reverse (const kmer_type& kmer)  { return revcomp (kmer, this->_sizeKmer); }

    /** */
    std::string toString (const kmer_type& kmer)  {  return kmer.toString (this->_sizeKmer);  }

    /** */
    std::pair<kmer_type, bool> getKmer (const tools::misc::Data& data, size_t idx=0, KmerMode mode = KMER_MINIMUM)
    {
        /** We compute the first kmer as a polynomial value. */
        kmer_type graine  = ModelAbstract<kmer_type>::codeSeed (data.getBuffer()+idx, data.getEncoding(), KMER_DIRECT);
        kmer_type revcomp = ModelAbstract<kmer_type>::codeSeed (data.getBuffer()+idx, data.getEncoding(), KMER_REVCOMP);

        bool      isForward = graine < revcomp;
        kmer_type value     = isForward ? graine : revcomp;

             if (mode == KMER_MINIMUM)  {  return std::make_pair (value, isForward); }
        else if (mode == KMER_DIRECT)   {  return std::make_pair (graine, true); }
        else  {  throw system::Exception ("BAD KMER MODE");  }
    }

    /** */
    bool build (tools::misc::Data& data, std::vector<kmer_type>& kmersBuffer,KmerMode mode = KMER_MINIMUM)
    {
        /** We compute the number of kmers for the provided data. Note that we have to check that we have
         * enough nucleotides according to the current kmer size. */
        int32_t nbKmers = data.size() - this->getSpan() + 1;
        if (nbKmers <= 0)  { return false; }

        /** We resize the resulting kmers vector. */
        kmersBuffer.resize (nbKmers);

        /** Shortcut used for computing kmers recursively. */
        char* buffer = data.getBuffer() + this->getSpan() - 1;

        /** We compute the first kmer as a polynomial value. */
        kmer_type graine  = ModelAbstract<kmer_type>::codeSeed (data.getBuffer(), data.getEncoding(), KMER_DIRECT);
        kmer_type revcomp = ModelAbstract<kmer_type>::codeSeed (data.getBuffer(), data.getEncoding(), KMER_REVCOMP);
        if(mode == KMER_MINIMUM)
        {
            kmersBuffer[0]    = std::min (graine, revcomp);
        }
        else if (mode == KMER_DIRECT)
        {
            kmersBuffer[0]    = graine;
        }
        else  {  throw system::Exception ("BAD KMER MODE");  }

        /** NOTE: we have some kind of duplicated code here: the only difference is the way we retrieve the 'c'
         *  character according toe the encoding mode. We could have some function pointer for factorizing this.
         *
         *  BUT, using function pointer would lead to many function calls with bad performance. We prefer therefore
         *  duplicate code and let inline work when needed.
         */
        if (data.getEncoding() == tools::misc::Data::ASCII && mode == KMER_MINIMUM)
        {
            /** We compute the next kmers in a recursive way. */
            for (size_t i=1; i<nbKmers; i++)
            {
                char c = ModelAbstract<kmer_type>::NT2int(buffer[i]);
                graine  = ( (graine  << 2) +  c                     ) & this->_kmerMask;
                revcomp = ( (revcomp >> 2) +  this->_revcompTable[c]) & this->_kmerMask;
                kmersBuffer[i] = std::min (graine, revcomp);
            }
        }
        else if (data.getEncoding() == tools::misc::Data::INTEGER && mode == KMER_MINIMUM)
        {
            /** We compute the next kmers in a recursive way. */
            for (size_t i=1; i<nbKmers; i++)
            {
                char c = ModelAbstract<kmer_type>::NTidentity(buffer[i]);
                graine  = ( (graine  << 2) +  c                     ) & this->_kmerMask;
                revcomp = ( (revcomp >> 2) +  this->_revcompTable[c]) & this->_kmerMask;
                kmersBuffer[i] = std::min (graine, revcomp);
            }
        }
        else if (data.getEncoding() == tools::misc::Data::BINARY && mode == KMER_MINIMUM)
        {
            buffer = data.getBuffer();

            size_t i0 = this->getSpan();
            size_t i1 = data.size();

            /** We compute the next kmers in a recursive way. */
            for (size_t i=i0; i<i1; i++)
            {
                char c = ((buffer[i>>2] >> ((3-(i&3))*2)) & 3);

                graine  = ( (graine  << 2) +  c                     ) & this->_kmerMask;
                revcomp = ( (revcomp >> 2) +  this->_revcompTable[c]) & this->_kmerMask;
                kmersBuffer[i-i0+1] = std::min (graine, revcomp);
            }
        }
        else if (data.getEncoding() == tools::misc::Data::ASCII && mode == KMER_DIRECT) //GR modif, check this with erwan
        {
            /** We compute the next kmers in a recursive way. */
            for (size_t i=1; i<nbKmers; i++)
            {
                char c = ModelAbstract<kmer_type>::NT2int(buffer[i]);
                graine  = ( (graine  << 2) +  c                     ) & this->_kmerMask;
                kmersBuffer[i] = graine;
            }
        }
        else if (data.getEncoding() == tools::misc::Data::INTEGER && mode == KMER_DIRECT)
        {
            /** We compute the next kmers in a recursive way. */
            for (size_t i=1; i<nbKmers; i++)
            {
                char c = ModelAbstract<kmer_type>::NTidentity(buffer[i]);
                graine  = ( (graine  << 2) +  c                     ) & this->_kmerMask;
                kmersBuffer[i] = graine;
            }
        }
        else if (data.getEncoding() == tools::misc::Data::BINARY && mode == KMER_DIRECT)
        {
            /** We compute the next kmers in a recursive way. */
            for (size_t i=1; i<nbKmers; i++)
            {
                char c = ((buffer[i>>2] >> ((3-(i&3))*2)) & 3);
                graine  = ( (graine  << 2) +  c                     ) & this->_kmerMask;
                kmersBuffer[i] = graine;
            }
        }
        else  {  throw system::Exception ("BAD DATA FORMAT TO BE CONVERTED IN KMERS...");  }

        return true;
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MODEL_HPP_ */
