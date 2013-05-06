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

/** \file Model.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Basic implementation of the IModel interface
 */

#ifndef _GATB_CORE_KMER_IMPL_MODEL_HPP_
#define _GATB_CORE_KMER_IMPL_MODEL_HPP_

/********************************************************************************/

#include <gatb/kmer/api/IModel.hpp>
#include <gatb/bank/api/IAlphabet.hpp>
#include <gatb/bank/impl/Alphabet.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/misc/api/Data.hpp>


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
template <typename kmer_type> class Model : public IModel<kmer_type>
{
public:

    /** Constructor.
     * \param[in] span : size of the kmers for this model
     * \param[in] alphabet : underlying alphabet
     */
    Model (size_t span, bank::IAlphabet& alphabet = bank::impl::AlphabetNucleic::singleton()) : _sizeKmer(span), _alphabet(alphabet)
    {
        _kmerMask = (((kmer_type)1) << (_sizeKmer*2))-1;
    }

    /** \copydoc IModel::getAlphabet */
    bank::IAlphabet&  getAlphabet() { return _alphabet; }

    /** \copydoc IModel::getSpan */
    size_t getSpan () { return _sizeKmer; }

    /** \copydoc IModel::getMemorySize */
    size_t getMemorySize ()  { return sizeof (kmer_type); }

    /** \copydoc IModel::codeSeed */
    kmer_type codeSeed (const char* seq, KmerMode mode)
    {
        kmer_type x = 0;
        for (size_t i=0; i<_sizeKmer; ++i)  {  x = x*4 + NT2int(seq[i]);  }
        //if (revcomp)  {  x = revcomp64(x); }
        return x;
    }

    /** \copydoc IModel::codeSeed */
    kmer_type codeSeed_bin (const char* seq, KmerMode mode)
    {
        kmer_type x = 0;
        for (size_t i=0; i<_sizeKmer; ++i)  {  x = x*4 + (seq[i]);  }
        //if (revcomp)  {  x = revcomp64(x); }
        return x;
    }

    /** \copydoc IModel::codeSeedRight
     * Copied from original function codeSeedRight from Minia */
    kmer_type codeSeedRight (const kmer_type& val_seed, char nucleotide, KmerMode mode)
    {
        switch (mode)
        {
        case KMER_REVCOMP:   return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[NT2int(nucleotide)]) <<  (2*(_sizeKmer-1))  )  ) & _kmerMask;
        case KMER_DIRECT:    return (val_seed*4 + NT2int (nucleotide) ) & _kmerMask;
        case KMER_BOTH:      return 0;
        default:        return 0;
        }
    }

    /** \copydoc IModel::codeSeedRight
     * Copied from original function codeSeedRight from Minia */
    kmer_type codeSeedRight_bin (const kmer_type& val_seed, char nucleotide, KmerMode mode)
    {
        switch (mode)
        {
        case KMER_REVCOMP:   return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[(nucleotide)]) <<  (2*(_sizeKmer-1))  )  ) & _kmerMask;
        case KMER_DIRECT:    return (val_seed*4 +  nucleotide) & _kmerMask;
        case KMER_BOTH:      return 0;
        default:        return 0;
        }
    }

    /************************************************************/
    /** \brief Specific Iterator impl for Model class */
    class Iterator : public tools::dp::Iterator<kmer_type>
    {
    public:
        /** Constructor.
         * \param[in] ref : the associated model instance.
         * \param[in] mode  : give the way kmers are computed
         */
        Iterator (Model& ref, KmerMode mode)
            : _ref(ref), _mode(mode), _idx(0), _idxMax(0), _isBinary(false), _bufferNext(0)    {  _data.buffer = _dataBuffer;   }

        /** Destructor */
        ~Iterator () {}

        /** Set the data to be iterated.
         * \param[in] data : the data as information source for the iterator
         */
        void setData (tools::misc::Data& data)
        {
            switch (data.getEncoding())
            {
            case tools::misc::Data::ASCII:    _data = data;  break;
            case tools::misc::Data::BINARY:
            {
                tools::misc::DataConverter::convert (data, _data);
                break;
            }
            default: break;
            }

            _idxMax     = _data.getSize()   - _ref._sizeKmer;
            _bufferNext = _data.getBuffer() + _ref._sizeKmer - 1;
            _isBinary   = _data.getEncoding() == tools::misc::Data::BINARY;
        }

        /** \copydoc tools::dp::Iterator::first */
        void first()
        {
            _idx = 0;
            if (isDone())  { return; }

            /** We compute the first kmer. The next kmer will be computed from this one. */
            if (_isBinary)   {  _current = _ref.codeSeed_bin (_data.getBuffer(), _mode);  }
            else             {  _current = _ref.codeSeed     (_data.getBuffer(), _mode);  }
        }

        /** \copydoc tools::dp::Iterator::next
         */
        void next()
        {
            if (_isBinary)  {  _current =  _ref.codeSeedRight_bin (_current, _bufferNext[++_idx], _mode);  }
            else            {  _current =  _ref.codeSeedRight     (_current, _bufferNext[++_idx], _mode);  }
        }

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _idx > _idxMax; }

        /** \copydoc tools::dp::Iterator::item */
        kmer_type& item ()     { return _current; }

    private:
        Model&      _ref;
        KmerMode    _mode;
        kmer_type   _current;
        tools::misc::Data _data;
        char        _dataBuffer[10*1024];
        size_t      _idx;
        size_t      _idxMax;
        bool        _isBinary;
        char*       _bufferNext;
    };

private:
    size_t           _sizeKmer;
    bank::IAlphabet& _alphabet;
    kmer_type        _kmerMask;

    /**  A=0, C=1, T=2, G=3 */
    int NT2int(char nt)
    {
        int i;
        i = nt;
        i = (i>>1)&3; // that's quite clever, guillaume.
        return i;
    }

    /** */
    u_int64_t revcomp64 (u_int64_t x)
    {
        u_int64_t revcomp = x;

        unsigned char * kmerrev  = (unsigned char *) (&revcomp);
        unsigned char * kmer     = (unsigned char *) (&x);

        for (size_t i=0; i<8; ++i)  {  kmerrev[7-i] = revcomp_4NT[kmer[i]];  }

        return (revcomp >> (2*( 4*sizeof(u_int64_t) - _sizeKmer))  ) ;
    }

    friend class Iterator;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MODEL_HPP_ */
