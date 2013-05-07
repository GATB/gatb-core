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

#include <gatb/system/api/Exception.hpp>
#include <gatb/kmer/api/IModel.hpp>
#include <gatb/bank/api/IAlphabet.hpp>
#include <gatb/bank/impl/Alphabet.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/misc/api/Data.hpp>

#include <vector>

#include <iostream>

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
    kmer_type codeSeed (const char* seq, KmerMode mode, int (*transfo)(char nucl) = &NT2int)
    {
        switch (mode)
        {
            case KMER_DIRECT:
            {
                kmer_type x = 0;
                for (size_t i=0; i<_sizeKmer; ++i)  {  x = x*4 + transfo(seq[i]);  }
                return x;
            }
            case KMER_REVCOMP:
            {
                // COULD BE OPTIMIZED
                kmer_type x = 0;
                for (size_t i=0; i<_sizeKmer; ++i)  {  x = x*4 + transfo(seq[i]);  }
                return revcomp64(x);
            }

            case KMER_MINIMUM:
            default:
            {
                // COULD BE OPTIMIZED
                kmer_type x = 0;
                for (size_t i=0; i<_sizeKmer; ++i)  {  x = x*4 + transfo(seq[i]);  }
                kmer_type y = revcomp64(x);
                return min (x,y);
            }
        }
    }

    /** \copydoc IModel::codeSeedRight
     * Copied from original function codeSeedRight from Minia */
    kmer_type codeSeedRight (const kmer_type& val_seed, char nucleotide, KmerMode mode, int (*transfo)(char nucl) = &NT2int)
    {
        switch (mode)
        {
        case KMER_DIRECT:    return (val_seed*4 + transfo (nucleotide) ) & _kmerMask;

        case KMER_REVCOMP:   return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[transfo(nucleotide)]) <<  (2*(_sizeKmer-1))  )  ) & _kmerMask;

        default:
        case KMER_MINIMUM:   return min (
                codeSeedRight (val_seed, nucleotide, KMER_DIRECT,  transfo),
                codeSeedRight (val_seed, nucleotide, KMER_REVCOMP, transfo)
            );
        }
    }

#if 1
    /************************************************************/
    /** \brief Specific Iterator impl for Model class */
    class Iterator : public tools::dp::Iterator<kmer_type>
    {
    public:
        /** Constructor.
         * \param[in] ref : the associated model instance.
         * \param[in] mode  : give the way kmers are computed
         */
        Iterator (Model& ref, KmerMode mode)  : _ref(ref), _mode(mode), _idx(0), _nbKmers(0)  {}

        /** Destructor */
        ~Iterator () {}

        /** Set the data to be iterated.
         *  Equivalent to the end of the Minia's KmersBuffer::readkmers()
         * \param[in] data : the data as information source for the iterator
         */
        void setData (tools::misc::Data& d)
        {
            /** By default, we will use the provided data with a ASCII encoding. */
            tools::misc::Data* data = &d;

            /** For ASCII encoding, we need to translate the nucleotides. */
            int (*transfo)(char nucl) = &NT2int;

            /** We may have to expand the binary data. */
            if (d.getEncoding() == tools::misc::Data::BINARY)
            {
                size_t expandedLen = d.getSize() ;
                if (_binaryVector.size () < expandedLen)
                {
                    _binaryVector.resize (expandedLen + 4);
                    _binaryData.buffer = _binaryVector.data();
                    _binaryData.size   = _binaryVector.size();
                }

                /** We convert the provided binary data into integer encoding. */
                tools::misc::DataConverter::convert (d, _binaryData);

                data = &_binaryData;

                /** For binary encoding, we need no nucleotide transformation. */
                transfo = &NTIdentity;
            }

            /** We compute the number of kmers for the provided data. */
            _nbKmers = data->getSize() - _ref.getSpan() + 1;

            /** We may have to resize the kmers buffer. Note the +1 => allow the 'next' method to go beyond
             * the real number of kmers, in order to escape one 'isDone' test. */
            if (_nbKmers >= _kmersBuffer.size())  { _kmersBuffer.resize (_nbKmers+1); }

            /** Shortcut used for computing kmers recursively. */
            char* buffer = data->getBuffer() + _ref.getSpan() - 1;

            if (_mode == KMER_DIRECT || _mode == KMER_REVCOMP)
            {
                /** We compute the first kmer as a polynomial value. */
                kmer_type graine = _ref.codeSeed (data->getBuffer(), _mode, transfo);
                _kmersBuffer[0] = graine;

                /** We compute the next kmers in a recursive way. */
                for (size_t i=1; i<_nbKmers; i++)
                {
                    graine = _ref.codeSeedRight  (graine, buffer[i], _mode, transfo);
                    _kmersBuffer[i] = graine;
                }
            }

            else if (_mode == KMER_MINIMUM)
            {
                /** We compute the first kmer as a polynomial value. */
                kmer_type graine  = _ref.codeSeed (data->getBuffer(), KMER_DIRECT,  transfo);
                kmer_type revcomp = _ref.codeSeed (data->getBuffer(), KMER_REVCOMP, transfo);
                _kmersBuffer[0] = min (graine, revcomp);

                /** We compute the next kmers in a recursive way. */
                for (size_t i=1; i<_nbKmers; i++)
                {
                    char c = buffer[i];

#if 1
                    graine  = _ref.codeSeedRight  (graine,  buffer[i], KMER_DIRECT,  transfo);
                    revcomp = _ref.codeSeedRight  (revcomp, buffer[i], KMER_REVCOMP, transfo);
#else
                    graine =  ( (graine << 2) +  c) & _ref._kmerMask;
                    revcomp = ((revcomp >> 2) +  ( ((kmer_type) comp_NT[c]) <<  (2*(_ref._sizeKmer-1))  )  ) & _ref._kmerMask;
#endif
                    _kmersBuffer[i] = min (graine, revcomp);
                }
            }
        }

        /** \copydoc tools::dp::Iterator::first */
        void first()
        {
            _idx = 0;
            _current = _kmersBuffer [_idx];
        }

        /** \copydoc tools::dp::Iterator::next */
        void next()  {  _current = _kmersBuffer [++_idx];  }

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  {  return _idx >= _nbKmers; }

        /** \copydoc tools::dp::Iterator::item */
        kmer_type& item ()     { return _current; }

    private:
        Model&      _ref;
        KmerMode    _mode;

        kmer_type   _current;
        std::vector<kmer_type> _kmersBuffer;

        int32_t _idx;
        int32_t _nbKmers;

        tools::misc::Data _binaryData;
        std::vector<char> _binaryVector;
    };

#else

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

#endif

private:
    size_t           _sizeKmer;
    bank::IAlphabet& _alphabet;
    kmer_type        _kmerMask;

    /**  A=0, C=1, T=2, G=3 */
    static int NT2int(char nt)
    {
        int i;
        i = nt;
        i = (i>>1)&3; // that's quite clever, guillaume.
        return i;
    }

    static int NTIdentity(char nt) { return nt; }

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
