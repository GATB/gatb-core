/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
template <typename kmer_type> class ModelAbstract : public IModel<kmer_type>, public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] span : size of the kmers for this model
     * \param[in] alphabet : underlying alphabet
     */
    ModelAbstract (size_t span, bank::IAlphabet& alphabet) : _sizeKmer(span), _alphabet(alphabet)
    {
    	/** We check that the kmer_type precision is enough for the required kmers span. */
    	if (span >= sizeof(kmer_type)*4)
    	{
    		throw system::Exception ("kmer_type '%s' has too low precision (%d bits) for the required %d kmer size",
				kmer_type().getName(), sizeof(kmer_type)*8, span
			);
    	}

    	/** We compute the mask of the kmer. Useful for computing kmers in a recursive way. */
    	kmer_type un = 1;
    	_kmerMask = (un << (_sizeKmer*2)) - un;

    	size_t shift = 2*(_sizeKmer-1);

        /** The _revcompTable is a shortcut used while computing revcomp recursively. */
        /** Important: don't forget the kmer_type cast, otherwise the result in only on 32 bits. */
        for (size_t i=0; i<4; i++)
        {
        	kmer_type tmp  = comp_NT[i];
        	_revcompTable[i] = tmp << shift;
        }
    }

    /** \copydoc IModel::getAlphabet */
    bank::IAlphabet&  getAlphabet() { return _alphabet; }

    /** \copydoc IModel::getSpan */
    size_t getSpan () { return _sizeKmer; }

    /** \return the kmer mask for the model, ie for the model kmer size. */
    kmer_type getMask () { return _kmerMask; }

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
         * \param[in] d : the data as information source for the iterator
         * \param[in] mode : mode for building the kmer
         */
        void setData (tools::misc::Data& d,KmerMode mode = KMER_MINIMUM)
        {
            /** We fill the vector with the items to be iterated. */
            if(mode == KMER_MINIMUM)
                _ref.build (d, this->_items);
            else
                _ref.build (d, this->_items,mode);

            /** We set the vector size. */
            this->_nb = this->_items.size();
        }

    private:
        /** Reference on the underlying model; called for its 'build' method. */
        ModelAbstract& _ref;
    };

    /************************************************************/
    class KmerNeighborIterator : public tools::dp::impl::VectorIterator<kmer_type>
    {
    public:

        /** */
        KmerNeighborIterator (ModelAbstract& ref)  : _ref(ref) {  this->_items.resize(8);   this->_nb = 8; }

        /** */
        virtual ~KmerNeighborIterator() {}

        /** */
        void setSource (const kmer_type& source)
        {
            size_t idx = 0;

            kmer_type rev = core::tools::math::revcomp (source, _ref.getSpan());

            /** We compute the 8 possible neighbors. */
            for (size_t nt=0; nt<4; nt++)
            {
                {
                    kmer_type next1 = (((source) * 4 )  + nt) & _ref.getMask();
                    kmer_type next2 = revcomp (next1, _ref.getSpan());
                    this->_items[idx++] = std::min (next1, next2);
                }
                {
                    kmer_type next1 = (((rev) * 4 )  + nt) & _ref.getMask();
                    kmer_type next2 = revcomp (next1, _ref.getSpan());
                    this->_items[idx++] = std::min (next1, next2);
                }
            }
        }

    private:
        /** Reference on the underlying model; called for its 'build' method. */
        ModelAbstract& _ref;
    };


    /** */
    template<typename Functor>
    void iterateNeighbors (const kmer_type& source, const Functor& fct)
    {
        kmer_type rev = core::tools::math::revcomp (source, getSpan());

        /** We compute the 8 possible neighbors. */
        for (size_t nt=0; nt<4; nt++)
        {
            {
                kmer_type next1 = (((source) * 4 )  + nt) & getMask();
                kmer_type next2 = revcomp (next1, getSpan());
                fct (std::min (next1, next2));
            }
            {
                kmer_type next1 = (((rev) * 4 )  + nt) & getMask();
                kmer_type next2 = revcomp (next1, getSpan());
                fct (std::min (next1, next2));
            }
        }
    }

protected:

    /** Size of a kmer for this model. */
    size_t           _sizeKmer;

    /** Alphabet telling what kind of letters are supposed to be used.
     * Not really used so far, but put here for the future if we will also deal with
     * amino acids alphabets for implementing other algorithms than assembly. */
    bank::IAlphabet& _alphabet;

    /** Mask for the kmer. Used for computing recursively kmers. */
    kmer_type  _kmerMask;

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

              if (encoding == tools::misc::Data::ASCII)    { for (size_t i=0; i<this->_sizeKmer; ++i)  {  direct = direct*4 + NT2int(seq[i]);  }  }
        else  if (encoding == tools::misc::Data::INTEGER)  { for (size_t i=0; i<this->_sizeKmer; ++i)  {  direct = direct*4 +       (seq[i]);  }  }
        else  if (encoding == tools::misc::Data::BINARY)   { for (size_t i=0; i<this->_sizeKmer; ++i)  {
                direct = direct*4 + (((seq[i>>2] >> ((3-(i&3))*2)) & 3));
            }
        }
        else  { throw system::Exception ("BAD FORMAT IN codeSeed"); }


        if (mode == KMER_DIRECT)  { return direct; }

        kmer_type rev = revcomp (direct, _sizeKmer);

        if (mode == KMER_REVCOMP)  { return rev; }

        return std::min (direct, rev);
    }

    /** Compute the successor of a kmer in a recursive way.
     *  WARNING ! we can't compute the minimum of direct and revcomp with this function since we need to know
     *  both current direct and revcomp for computing the next one.
     * \param[in] seed : initial kmer from which we want to compute the successor
     * \param[in] nucl : the nucleotide to be appended to the current kmer
     * \param[in] encoding : encoding of the source data
     * \param[in] mode : tells how to compute the kmer.
     */
    kmer_type codeSeedRight (const kmer_type& seed, char nucl, tools::misc::Data::Encoding_e encoding, KmerMode mode)
    {
        kmer_type direct;

        if (encoding == tools::misc::Data::ASCII)   {  direct = ( (seed << 2) + NT2int(nucl)) & _kmerMask;  }
        else                                        {  direct = ( (seed << 2) +       (nucl)) & _kmerMask;  }

        if (mode == KMER_DIRECT)  { return direct; }

        kmer_type rev = revcomp (direct, _sizeKmer);

        if (mode == KMER_REVCOMP)  { return rev; }

        if (mode == KMER_MINIMUM)  { return std::min(direct,rev); }
        
        throw system::Exception ("BAD MODE for computing kmer successor");
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MODEL_ABSTRACT_HPP_ */
