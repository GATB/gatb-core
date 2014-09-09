/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/** \file Model.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Kmer management
 */

#ifndef _GATB_CORE_KMER_IMPL_MODEL_HPP_
#define _GATB_CORE_KMER_IMPL_MODEL_HPP_

/********************************************************************************/

#include <gatb/system/api/Exception.hpp>
#include <gatb/kmer/api/IModel.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/Data.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>

#include <gatb/tools/math/LargeInt.hpp>

#include <vector>
#include <algorithm>
#include <iostream>

extern const char bin2NT[] ;
extern const char binrev[] ;
extern const unsigned char revcomp_4NT[];
extern const unsigned char comp_NT    [];

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

#define KMER_DEFAULT_SPAN KSIZE_1

/********************************************************************************/

/** \brief Entry point for kmer management.
 *
 * This structure is only a container for other types defined inside. The specificity is
 * that this structure is templated by a 'span' integer that represents the maximal kmer
 * size supported (actually, the max value is 'span-1').
 *
 * Inside this structure, we have the following main elements:
 *      - 'Type'  : this is the integer type representing kmer values
 *      - 'Model' : provides many services for managing kmers
 *      - 'Count' : roughly speaking, this a kmer value with an associated abundance
 *
 * This structure must be used only with for 4 values (32,64,96,128 for instance, see Model.cpp), otherwise
 * a compilation error occurs (more values could be added in the future).
 *
 * A default value of 32 is defined for the template parameter, so writing 'Kmer<>::Model'
 * represents a model that supports kmers of size up to 31 (included).
 */
template <size_t span=KMER_DEFAULT_SPAN>
struct Kmer
{
    /************************************************************/
    /***********************     TYPE     ***********************/
    /************************************************************/

    /** Alias type for the integer value of a kmer. We use the LargeInt class for supporting big integers.
     * Note that the template parameter 'span' represents the maximal kmer size supported by the Kmer class.
     * A conversion to the template parameter of LargeInt is done.
     */
    typedef tools::math::LargeInt<(span+31)/32> Type;


    /************************************************************/
    /***********************     MODEL    ***********************/
    /************************************************************/

    /** Forward declarations. */
    class ModelDirect;
    class ModelCanonical;
    template<class Model, class Comparator> class ModelMinimizer;

    /** Now, we need to define what is a kmer for each kind of model.
     *
     * The simple case is KmerDirect, where only the value of the kmer is available
     * as a method 'value' returning a Type object.
     *
     * The second case is KmerCanonical, which is the same as KmerDirect, but with
     * two other methods 'forward' and 'revcomp'
     *
     * The third case is KmerMinimizer<Model> which allows to handle minimizers associated
     * to a kmer. This class inherits from the Model::Kmer type and adds methods specific
     * to minimizers, such as 'minimizer' itself (ie the Model::Kmer object holding the
     * minimizer), 'position' giving the position of the minimizer whithin the kmer and
     * 'hasChanged' telling whether a minimizer has changed during iteration of kmers from
     * some data source (a sequence data for instance).
     */

    /** Kmer type for the ModelDirect class. */
    class KmerDirect
    {
    public:
        /** Returns the value of the kmer.
         * \return the kmer value as a Type object. */
        const Type& value  () const { return _value;   }

        /** Comparison operator between two instances.
         * \param[in] t : object to be compared to
         * \return true if the values are the same, false otherwise. */
        bool operator< (const KmerDirect& t) const  { return this->_value < t._value; };

        /** Set the value of the kmer
         * \param[in] val : value to be set. */
        void set (const Type& val) { _value=val; }

    protected:
        Type _value;
        friend class ModelDirect;
    };

    /** Kmer type for the ModelCanonical class. */
    class KmerCanonical : public KmerDirect
    {
    public:

        /** Set the value attribute. */
        void set (const Type& value)  {  KmerDirect::set(value); }

        /** Set the forward/revcomp attributes. */
        void set (const Type& forward, const Type& revcomp)
        {
            _forward=forward;
            _revcomp=revcomp;
            KmerDirect::set (std::min (_forward,_revcomp));
        }

        /** Returns the forward value of this canonical kmer.
         * \return the forward value */
        const Type& forward() const { return _forward; }

        /** Returns the reverse complement value of this canonical kmer.
         * \return the reverse complement value */
        const Type& revcomp() const { return _revcomp; }

    protected:
        Type _forward;  Type _revcomp;
        friend class ModelCanonical;
    };

    /** Kmer type for the ModelMinimizer class. */
    template<class Model, class Comparator>
    class KmerMinimizer : public Model::Kmer
    {
    public:

        /** Returns the minimizer of the current kmer as a Model::Kmer object
         * \return the Model::Kmer instance */
        const typename Model::Kmer& minimizer() const  {  return minimizers[minimizerIdx];  }

        /** Returns the position of the minimizer within the kmer.
         * \return the position of the minimizer. */
        size_t position () const
        {
            /** By convention, if there is no minimizer, we return position 0. */
            if (isDefined()==false) { return 0; }

            return startIdx<minimizerIdx ? minimizerIdx-startIdx-1 : (minimizerIdx+nbMinimizer)-startIdx-1;
        }

        /** Tells whether the minimizer has changed; useful while iterating kmers
         * \return true if changed, false otherwise */
        bool hasChanged () const  {  return changed;  }

        /** Tells whether the minimizer is defined within the kmer.
         * \return true if defined, false otherwise. */
        bool isDefined () const { return minimizerIdx!=nbMinimizer; }

    protected:


        typename Model::Kmer minimizers[span];
        size_t minimizerIdx;
        size_t startIdx;
        size_t nbMinimizer;
        bool changed;
        friend class ModelMinimizer<Model,Comparator>;
    };

    /** Abstract class that provides kmer management.
     *
     * This class is the base class for kmer management. It provides several services on this purpose
     * like getting kmer information from some nucleotides sequence, or iterate kmers through such
     * a sequence.
     *
     * This class has two templates types :
     *
     *      1) ModelImpl : ModelAbstract is design for static polymorphism and ModelImpl is the implementation
     *                     that must be provided to it
     *
     *      2) T : type of kmers handled by the class (ie KmerDirect, KmerCanonical...); I was not successful
     *             in trying to hide KmerXXX classes in the dedicated ModelXXX classes because of mutual
     *             dependencies while template specializations (maybe a solution one day)
     *
     * End user will be given instances of Kmer class, delivering more or less information according to the
     * specific type of ModelImpl
     */
    template <class ModelImpl, typename T>
    class ModelAbstract : public system::SmartPointer
    {
    public:

        /** Type of kmers provided by the class. */
        typedef T Kmer;

        /** (default) Constructor. The provided (runtime) kmer size must be coherent with the span (static) value.
         * \param[in] sizeKmer : size of kmers handled by the instance.*/
        ModelAbstract (size_t sizeKmer=span-1) : _kmerSize(sizeKmer)
        {
            /** We check that the Type precision is enough for the required kmers span. */
            if (sizeKmer >= span)
            {
                throw system::Exception ("Type '%s' has too low precision (%d bits) for the required %d kmer size",
                    Type().getName(), Type().getSize(), sizeKmer
                );
            }

            /** We compute the mask of the kmer. Useful for computing kmers in a recursive way. */
            Type un = 1;
            _kmerMask = (un << (_kmerSize*2)) - un;

            size_t shift = 2*(_kmerSize-1);

            /** The _revcompTable is a shortcut used while computing revcomp recursively. */
            /** Important: don't forget the Type cast, otherwise the result in only on 32 bits. */
            for (size_t i=0; i<4; i++)   {  Type tmp  = comp_NT[i];  _revcompTable[i] = tmp << shift;  }
        }

        /** Returns the span of the model
         * \return the model span. */
        size_t getSpan () const { return span; }

        /** Get the memory size (in bytes) of a Kmer<span>::Type object.
         * \return the memory size of a kmer. */
        size_t getMemorySize ()  const  { return sizeof (Type); }

        /** Gives the kmer size for this model.
         * \return the kmer size. */
        size_t getKmerSize () const { return _kmerSize; }

        /** Gives the maximum value of a kmer for the instance.
         * \return the maximum kmer value. */
        const Type& getKmerMax () const { return _kmerMask; }

        /** Returns an ascii representation of the kmer value.
         * \param[in] kmer : the kmer we want an ascii representation for
         * \return a string instance holding the ascii representation. */
        std::string toString (const Type& kmer) const  {  return kmer.toString(_kmerSize);  }

        /** Compute the reverse complement of a kmer.
         * \param[in] kmer : the kmer to be reverse-completed.
         * \return the reverse complement. */
        Type reverse (const Type& kmer)  const  { return revcomp (kmer, this->_kmerSize); }

        /** Build a kmer from a Data object (ie a sequence of nucleotides), starting at an index in the nucleotides sequence.
         * The result is a pair holding the built kmer and a boolean set to yes if the built kmer has to be understood in
         * the forward sense, false otherwise.
         * \param[in] data : the data from which we extract a kmer
         * \param[in] idx : start index in the data object (default to 0)
         * \return a pair with the built kmer and a boolean set to yes if the kmer is understood in the forward strand
         */
        Kmer getKmer (const tools::misc::Data& data, size_t idx=0)  const
        {
            return codeSeed (data.getBuffer() + idx, data.getEncoding());
        }

        /** Iteration of the kmers from a data object through a functor (so lambda expressions can be used).
         * \param[in] data : the sequence of nucleotides.
         * \param[in] fct  : functor that handles one kmer */
        template<typename Callback>
        bool iterate (tools::misc::Data& data, Callback callback) const
        {
            return execute <Functor_iterate<Callback> > (data.getEncoding(), Functor_iterate<Callback>(data,callback));
        }

        /** Compute the kmer given some nucleotide data.
         *  Note that we don't check if we have enough nucleotides in the provided data.
         * \param[in] seq : the sequence
         * \param[in] encoding : encoding mode of the sequence
         * \return the kmer for the given nucleotides. */
        Kmer codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding) const
        {
            return execute<Functor_codeSeed> (encoding, Functor_codeSeed(seq));
        }

        /** Compute the next right kmer given a current kmer and a nucleotide.
         * \param[in] kmer : the current kmer as a starting point
         * \param[in] nucl : the next nucleotide
         * \param[in] encoding : encoding mode of the sequence
         * \return the kmer on the right of the given kmer. */
        Kmer codeSeedRight (const Kmer& kmer, char nucl, tools::misc::Data::Encoding_e encoding)  const
        {
            return execute<Functor_codeSeedRight> (encoding, Functor_codeSeedRight(kmer,nucl));
        }

        /** Build a vector of successive kmers from a given sequence of nucleotides provided as a Data object.
         * \param[in] data : the sequence of nucleotides.
         * \param[out] kmersBuffer : the successive kmers built from the data object.
         * \return true if kmers have been extracted, false otherwise. */
        bool build (tools::misc::Data& data, std::vector<Kmer>& kmersBuffer)  const
        {
            /** We compute the number of kmers for the provided data. Note that we have to check that we have
             * enough nucleotides according to the current kmer size. */
            int32_t nbKmers = data.size() - this->getKmerSize() + 1;
            if (nbKmers <= 0)  { return false; }

            /** We resize the resulting kmers vector. */
            kmersBuffer.resize (nbKmers);

            /** We fill the vector through a functor. */
            this->iterate (data, BuildFunctor<Kmer>(kmersBuffer));

            return true;
        }

        /** Iterate the neighbors of a given kmer; these neighbors are:
         *  - 4 outcoming neighbors
         *  - 4 incoming neighbors.
         *  This method uses a functor that will be called for each possible neighbor of the source kmer.
         *  \param[in] source : the kmer from which we want neighbors.
         *  \param[in] fct : a functor called for each neighbor.*/
        template<typename Functor>
        void iterateNeighbors (const Type& source, const Functor& fct)  const
        {
            Type rev = core::tools::math::revcomp (source, getKmerSize());

            /** We compute the 8 possible neighbors. */
            for (size_t nt=0; nt<4; nt++)
            {
                {
                    Type next1 = (((source) * 4 )  + nt) & getKmerMax();
                    Type next2 = revcomp (next1, getKmerSize());
                    fct (std::min (next1, next2));
                }
                {
                    Type next1 = (((rev) * 4 )  + nt) & getKmerMax();
                    Type next2 = revcomp (next1, getKmerSize());
                    fct (std::min (next1, next2));
                }
            }
        }

        /************************************************************/
        /** \brief Iterator on successive kmers
         *
         * This class will iterate successive kmers extracted from a Data object.
         * It is similar to the Model::build, except that here we don't have a container
         * holding all the successive kmers (ie. we have here only sequential access and
         * not direct access).
         *
         * To be used, such an iterator must be initialized with some sequence of nucleotides,
         * which is done with the 'setData' method.
         */
        class Iterator : public tools::dp::impl::VectorIterator<Kmer>
        {
        public:
            /** Constructor.
             * \param[in] ref : the associated model instance.
             */
            Iterator (ModelAbstract& ref)  : _ref(ref)   {}

            /** Set the data to be iterated.
             * \param[in] d : the data as information source for the iterator
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
        size_t  _kmerSize;

        /** Mask for the kmer. Used for computing recursively kmers. */
        Type  _kmerMask;

        /** Shortcut for easing/speeding up the recursive revcomp computation. */
        Type _revcompTable[4];

        /** */
        struct ConvertASCII    { static char get (const char* buffer, size_t idx)  { return (buffer[idx]>>1) & 3; }};
        struct ConvertInteger  { static char get (const char* buffer, size_t idx)  { return buffer[idx]; }         };
        struct ConvertBinary   { static char get (const char* buffer, size_t idx)  { return ((buffer[idx>>2] >> ((3-(idx&3))*2)) & 3); } };

        /** */
        template<class Convert>
        void polynom (const char* seq, Type& kmer)  const
        {  kmer = 0;  for (size_t i=0; i<_kmerSize; ++i)  {  kmer = (kmer<<2) + Convert::get(seq,i);  }  }

        /** Generic function that switches to the correct implementation according to the encoding scheme
         * of the provided Data parameter; the provided functor class is specialized with the correct data conversion type
         * and the called.
         */
        template<class Functor>
        typename Functor::Result execute (tools::misc::Data::Encoding_e encoding, Functor action) const
        {
            switch (encoding)
            {
                case  tools::misc::Data::ASCII:    return action.template operator()<ConvertASCII  > (this);
                case  tools::misc::Data::INTEGER:  return action.template operator()<ConvertInteger> (this);
                case  tools::misc::Data::BINARY:   return action.template operator()<ConvertBinary>  (this);
                default:  throw system::Exception ("BAD FORMAT IN 'execute'");
            }
        }

        /** Adaptor between the 'execute' method and the 'codeSeed' method. */
        struct Functor_codeSeed
        {
            typedef typename ModelImpl::Kmer Result;
            const char* buffer;
            Functor_codeSeed (const char* buffer) : buffer(buffer) {}
            template<class Convert>  Result operator() (const ModelAbstract* model)
            {
                Result result;
                static_cast<const ModelImpl*>(model)->template first <Convert> (buffer, result);
                return result;
            }
        };

        /** Adaptor between the 'execute' method and the 'codeSeedRight' method. */
        struct Functor_codeSeedRight
        {
            typedef typename ModelImpl::Kmer Result;
            const Kmer& kmer; char nucl;
            Functor_codeSeedRight (const Kmer& kmer, char nucl) : kmer(kmer), nucl(nucl) {}
            template<class Convert>  Result operator() (const ModelAbstract* model)
            {
                Result result=kmer;
                static_cast<const ModelImpl*>(model)->template next <Convert> (Convert::get(&nucl,0), result);
                return result;
            }
        };

        /** Adaptor between the 'execute' method and the 'iterate' method. */
        template<class Callback>
        struct Functor_iterate
        {
            typedef bool Result;
            tools::misc::Data& data; Callback callback;
            Functor_iterate (tools::misc::Data& data, Callback callback) : data(data), callback(callback) {}
            template<class Convert>  Result operator() (const ModelAbstract* model)
            {
                return static_cast<const ModelImpl*>(model)->template iterate<Callback, Convert> (data.getBuffer(), data.size(), callback);
            }
        };

        /** Template method that iterates the kmer of a given Data instance.
         *  Note : we use static polymorphism here (http://en.wikipedia.org/wiki/Template_metaprogramming)
         */
        template<typename Callback, typename Convert>
        bool iterate (const char* seq, size_t length, Callback callback) const
        {
            /** We compute the number of kmers for the provided data. Note that we have to check that we have
             * enough nucleotides according to the current kmer size. */
            int32_t nbKmers = length - _kmerSize + 1;
            if (nbKmers <= 0)  { return false; }

            /** We create a result instance. */
            typename ModelImpl::Kmer result;

            /** We compute the initial seed from the provided buffer. */
            static_cast<const ModelImpl*>(this)->template first<Convert> (seq, result);

            /** We need to keep track of the computed kmers. */
            size_t idxComputed = 0;

            /** We notify the result. */
            this->notification<Callback> (result, idxComputed, callback);

            /** We compute the following kmers from the first one.
             * We have consumed 'kmerSize' nucleotides so far for computing the first kmer,
             * so we start the loop with idx=_kmerSize.
             */
            for (size_t idx=_kmerSize; idx<length; idx++)
            {
                /** We get the current nucleotide. */
                char c = Convert::get (seq, idx);

                /** We compute the next kmer from the previous one. */
                static_cast<const ModelImpl*>(this)->template next<Convert> (c, result);

                /** We notify the result. */
                this->notification<Callback> (result, ++idxComputed, callback);
            }

            return true;
        }

        template <class Callcack>
        void  notification (const Kmer& value, size_t idx, Callcack callback) const {  callback (value, idx);  }

        /** */
        template<typename Type>
        struct BuildFunctor
        {
            std::vector<Type>& kmersBuffer;
            BuildFunctor (std::vector<Type>& kmersBuffer) : kmersBuffer(kmersBuffer) {}
            void operator() (const Type& kmer, size_t idx)  {  kmersBuffer[idx] = kmer;  }
        };
    };

    /********************************************************************************/

    /** \brief Model that handles "direct" kmers, ie sequences of nucleotides.
     * The associated value of such a kmer is computed as a polynom P(X) with X=4
     * and where the coefficients are in [0..3].
     * By convention, we use A=0, C=1, T=2 and G=3
     */
    class ModelDirect :  public ModelAbstract<ModelDirect, Kmer<span>::KmerDirect>
    {
    public:

        /** Type holding all the information of a kmer.  */
        typedef Kmer<span>::KmerDirect Kmer;

        /** Constructor.
         * \param[in] kmerSize : size of the kmers handled by the model. */
        ModelDirect (size_t kmerSize=span-1) : ModelAbstract<ModelDirect, Kmer> (kmerSize) {}

        /** Computes a kmer from a buffer holding nucleotides encoded in some format.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] buffer : holds the nucleotides sequence from which the kmer has to be computed
         * \param[out] value : kmer as a result
         */
        template <class Convert>
        void first (const char* buffer, Kmer& value)   const
        {
            this->template polynom<Convert> (buffer, value._value);
        }

        /** Computes a kmer in a recursive way, ie. from a kmer and the next
         * nucleotide. The next nucleotide is computed from a buffer and the
         * index of the nucleotide within the buffer.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] buffer : holds the nucleotides sequence from which the kmer has to be computed
         * \param[in] idx : index of the nucleotide within the buffer
         * \param[out] value kmer as a result
         */
        template <class Convert>
        void  next (char c, Kmer& value)   const
        {
            value._value = ( (value._value << 2) +  c) & this->_kmerMask;
        }
    };

    /********************************************************************************/

    /** \brief Model that handles "canonical" kmers, ie the minimum value of the
     * direct kmer and its reverse complement.
     */
    class ModelCanonical :  public ModelAbstract<ModelCanonical, Kmer<span>::KmerCanonical>
    {
    public:

        /** Type holding all the information of a kmer.  */
        typedef Kmer<span>::KmerCanonical Kmer;

        /** Constructor.
         * \param[in] kmerSize : size of the kmers handled by the model. */
        ModelCanonical (size_t kmerSize=span-1) : ModelAbstract<ModelCanonical, Kmer> (kmerSize) {}

        /** Computes a kmer from a buffer holding nucleotides encoded in some format.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] buffer : holds the nucleotides sequence from which the kmer has to be computed
         * \param[out] value : kmer as a result
         */
        template <class Convert>
        void  first (const char* seq, Kmer& value)   const
        {
            this->template polynom<Convert> (seq, value._forward);
            value._revcomp = this->reverse (value._forward);
            value._value   = std::min (value._forward, value._revcomp);
        }

        /** Computes a kmer in a recursive way, ie. from a kmer and the next
         * nucleotide. The next nucleotide is computed from a buffer and the
         * index of the nucleotide within the buffer.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] buffer : holds the nucleotides sequence from which the kmer has to be computed
         * \param[in] idx : index of the nucleotide within the buffer
         * \param[out] value kmer as a result
         */
        template <class Convert>
        void  next (char c, Kmer& value)   const
        {
            value._forward = ( (value._forward << 2) +  c) & this->_kmerMask;
            value._revcomp = ( (value._revcomp >> 2) +  this->_revcompTable[c]) & this->_kmerMask;
            value._value   = std::min (value._forward, value._revcomp);
        }
    };

    /********************************************************************************/

    struct ComparatorMinimizer
    {
        template<class Model>  void init (const Model& model, Type& best) const { best = model.getKmerMax(); }
        bool operator() (const Type& current, const Type& best) const { return current < best; }
    };

    /** \brief Model that handles kmers of the Model type + a minimizer
     *
     * This model supports the concept of minimizer. It acts as a Model instance (given as a
     * template class) and add minimizer information to the Kmer type.
     */
    template<class ModelType, class Comparator=ComparatorMinimizer>
    class ModelMinimizer :  public ModelAbstract <ModelMinimizer<ModelType,Comparator>, KmerMinimizer<ModelType,Comparator> >
    {
    public:

        /** Type holding all the information of a kmer.  */
        typedef KmerMinimizer<ModelType,Comparator> Kmer;

        /** Return a reference on the model used for managing mmers. */
        const ModelType& getMmersModel() const { return _miniModel; }

        /** Constructor.
         * \param[in] kmerSize      : size of the kmers handled by the model.
         * \param[in] minimizerSize : size of the mmers handled by the model. */
        ModelMinimizer (size_t kmerSize, size_t minimizerSize, Comparator cmp=Comparator())
            : ModelAbstract <ModelMinimizer<ModelType,Comparator>, Kmer > (kmerSize),
              _kmerModel(kmerSize), _miniModel(minimizerSize), _nbMinimizers(0), _cmp(cmp)
        {
            if (kmerSize <= minimizerSize)  { throw system::Exception ("Bad values for kmer %d and minimizer %d", kmerSize, minimizerSize); }

            /** We compute the number of mmers found in a kmer. */
            _nbMinimizers = _kmerModel.getKmerSize() - _miniModel.getKmerSize() + 1;
        }

        template <class Convert>
        void  first (const char* seq, Kmer& value)   const
        {
            /** We initialize a sentinel value at the end of the mmers vector.
             * This value will be used in case no minimizer is found.
             * The value is actually set by the Comparator instance provided as a template of the class. */
            _cmp.template init<ModelType> (getMmersModel(), (Type&)(value.minimizers [_nbMinimizers].value()));

            /** We memorize the number of minimizers. */
            value.nbMinimizer = _nbMinimizers;

            /** We compute the kmer. */
            _kmerModel.template first<Convert> (seq, value);

            /** We compute N potential minimizers and put them into the circular buffer. */
            _miniModel.template first<Convert> (seq, value.minimizers[0]);

            /** We have to shift the current buffer to extract next mmers. */
            seq += _miniModel.getKmerSize() - 1;

            for (size_t idx=1; idx<_nbMinimizers; idx++)
            {
                value.minimizers[idx] = value.minimizers[idx-1];
                _miniModel.template next<Convert> (Convert::get (seq, idx), value.minimizers[idx]);
            }

            /** We initialize the circular buffer index. */
            value.startIdx = _nbMinimizers - 1;

            /** We get the index of the the minimizer in the circular buffer. */
            value.minimizerIdx = getMinimizerIdx (value);
            value.changed   = true;
        }

        template <class Convert>
        void  next (char c, Kmer& value)   const
        {
            /** We get a copy of the current minimizer. */
            Type currentMinimizer = value.minimizers[value.minimizerIdx].value();

            size_t nextIdx =  (value.startIdx + 1) % _nbMinimizers;
            value.minimizers[nextIdx] = value.minimizers[value.startIdx];

            _kmerModel.template next<Convert> (c, value);
            _miniModel.template next<Convert> (c, value.minimizers[nextIdx]);

            /** By default, we set the minimizer has unchanged for the 'next' call. */
            value.changed = false;

            /** We update the starting index in the circular buffer. */
            value.startIdx = nextIdx;

            /** We may have to update the minimizer index in the following cases :
             *      1) the minimizer index is invalid
             *      2) the minimizer index is out of the current kmer window
             *      3) the new current mmer is best than the current minimizer
             */
            if (value.minimizerIdx==_nbMinimizers || value.minimizerIdx==nextIdx
                || _cmp (value.minimizers[nextIdx].value(), value.minimizers[value.minimizerIdx].value()) )
            {
                /** We update the minimizer index. */
                value.minimizerIdx = getMinimizerIdx (value);

                /** We check whether it is a true minimizer change. Note: checking only the indexes is not enough
                 * because we can have the same minimizer twice or more in the same kmer. */
                if (currentMinimizer != value.minimizers[value.minimizerIdx].value())  {  value.changed = true;  }
            }
        }

    private:
        ModelType _kmerModel;
        ModelType _miniModel;

        size_t _nbMinimizers;

        Comparator _cmp;

        /** Returns the minimizer of the provided vector of mmers. */
        int getMinimizerIdx(const Kmer& kmer) const
        {
            int result = _nbMinimizers;
            typename ModelType::Kmer current;  current.set(kmer.minimizers[result].value());

            /** we have to loop nbMinimizers but not starting from startIdx instead of 0
             *  => we split the loop in two parts. */
            size_t i0 = kmer.startIdx + 1;

            for (size_t i=i0; i<_nbMinimizers; i++)
            {
                if (_cmp(kmer.minimizers[i].value(), current.value())==true)  {  current = kmer.minimizers [result = i];  }
            }

            for (size_t i=0; i<i0; i++)
            {
                if (_cmp(kmer.minimizers[i].value(), current.value())==true)  {  current = kmer.minimizers [result = i];  }
            }

            /** We return the result. */
            return result;
        }
    };

    /************************************************************/
    /***********************     COUNT    ***********************/
    /************************************************************/
    /** \brief Structure associating a kmer value with an abundance value.
     *
     * This structure is useful for methods that counts kmer, such as the SortingCount algorithm.
     * It is also interesting to save [kmer,abundance] in a HDF5 format.
     *
     * By default, the abundance value is coded on 16 bits, so abundance up to 1<<16 can be used.
     */
    struct Count : tools::misc::Abundance<Type,u_int16_t>
    {
        /** Shortcut. */
        typedef u_int16_t Number;

        /** Constructor.
         * \param[in] val : integer value of the kmer
         * \param[in] abund : abundance for the kmer */
        Count(const Type& val, const Number& abund) : tools::misc::Abundance<Type,Number>(val, abund) {}

        /** Default constructor. */
        Count() : tools::misc::Abundance<Type,Number>(Type(), 0) {}

        /** Copy constructor. */
        Count(const Count& val) : tools::misc::Abundance<Type,Number>(val.value, val.abundance) {}

        /** Comparison operator
         * \param[in] other : object to be compared to
         * \return true if the provided kmer value is greater than the current one. */
        bool operator< (const Count& other) const {  return this->value < other.value; }
        
        /** Equal operator
         * \param[in] other : object to be compared to
         * \return true if the provided kmer value is greater than the current one. */
        bool operator== (const Count& other) const {  return (this->value == other.value && this->abundance == other.abundance); }
    };

};  // struct Kmer

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MODEL_HPP_ */
