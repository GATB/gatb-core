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

#define KMER_DEFAULT_SPAN 32

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
 * This structure must be used only with for 4 values (32,64,96,128, see Model.cpp), otherwise
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

    /** \brief Class providing services for kmer management.
     *
     * The Model class provides means to generate kmers, for instance from a sequence of nucleotides.
     *
     * It is defined by the kmer size (provided at construction). If the provided kmer size is not compatible
     * with the 'span' template parameter, an exception is thrown.
     *
     * The way kmers are generated through the Model class is defined by a KmerMode; an integer value of
     * type 'Type' is associated to a kmer, an this integer value may be computed differently according to
     * the chosen KmerMode.
     */
    class Model : public system::SmartPointer
    {
    public:

        /** Constructor.
         * \param[in] sizeKmer : size of the kmers for this model. If not set, use the maximum size allowed by the template parameter
         * \param[in] mode : tell how to compute kmer; default mode is 'minimum', ie min(forward,revcomp)
         */
        Model (size_t sizeKmer=span-1, KmerMode mode=KMER_MINIMUM);

        /** Gives the span, ie the maximum size of kmers allowed by this model. For instance, for span=32, the kmer size
         * can be as large as 31 (included). Actually, only a few values of span are allowed (32, 64, 96 and 128), otherwise
         * a compilation error will happen.
         * \return the span of this model.
         * */
        size_t getSpan () const { return span; }

        /** Gives the kmer size for this model.
         * \return the kmer size. */
        size_t getKmerSize () const { return _kmerSize; }

        /** Returns the maximum kmer value allowed for this model according to the provided kmer size.
         * For instance, if kmersize=31, the mask will be (1 << 2*kmersize) - 1
         * \return the kmer mask for the model, ie for the model kmer size. */
        Type getMask ()  const { return _kmerMask; }

        /**  Returns the kmer mode for this model.
         * \return the kmer mask for the model, ie for the model kmer size. */
        KmerMode getMode ()  const  { return _mode; }

        /** Get the memory size (in bytes) of a Kmer<span>::Type object.
         * \return the memory size of a kmer. */
        size_t getMemorySize ()  const  { return sizeof (Type); }

        /** Get a human readable string for the kmer (provided as an integer of type Type).
         * \return the string as a sequence of nucleotides. */
        std::string toString (const Type& kmer) const  {  return kmer.toString(_kmerSize);  }

        /** Compute the kmer given some nucleotide data.
         *  Note that we don't check if we have enough nucleotides in the provided data.
         * \param[in] seq : the sequence
         * \param[in] encoding : encoding mode of the sequence
         * \return the kmer for the given nucleotides. */
        Type codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding)   const  {  return codeSeed (seq, encoding, _mode);  }

        /** Compute the next right kmer given a current kmer and a nucleotide.
         * \param[in] seed : the current kmer as a starting point
         * \param[in] nucl : the next nucleotide
         * \param[in] encoding : encoding mode of the sequence
         * \return the kmer on the right of the given kmer. */
        Type codeSeedRight (const Type& seed, char nucl, tools::misc::Data::Encoding_e encoding)  const  {  return codeSeedRight (seed, nucl, encoding, _mode);  }

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
        std::pair<Type, bool> getKmer (const tools::misc::Data& data, size_t idx=0)  const  {  return getKmer (data, idx, _mode);  }

        /** Build a vector of successive kmers from a given sequence of nucleotides provided as a Data object.
         * \param[in] data : the sequence of nucleotides.
         * \param[out] kmersBuffer : the successive kmers built from the data object.
         * \return true if kmers have been extracted, false otherwise. */
        bool build (tools::misc::Data& data, std::vector<Type>& kmersBuffer)  const  {  return build (data, kmersBuffer, _mode);  }

        /** Iteration of the kmers from a data object through a functor (so lambda expressions can be used).
         * \param[in] data : the sequence of nucleotides.
         * \param[in] fct  : functor that handles one kmer */
        template<typename Functor> bool iterate (tools::misc::Data& data, Functor fct) const { return iterate (data, fct, _mode); }

        /** Iteration of all possible kmers through a functor (so lambda expressions can be used).
         * \param[in] fct  : functor that handles one kmer */
        template<typename Functor> bool iterate (Functor fct) const;

        /** Iterate the neighbors of a given kmer; these neighbors are:
         *  - 4 outcoming neighbors
         *  - 4 incoming neighbors.
         *  This method uses a functor that will be called for each possible neighbor of the source kmer.
         *  \param[in] source : the kmer from which we want neighbors.
         *  \param[in] fct : a functor called for each neighbor.*/
        template<typename Functor>
        void iterateNeighbors (const Type& source, const Functor& fct)  const;

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
        class Iterator : public tools::dp::impl::VectorIterator<Type>
        {
        public:
            /** Constructor.
             * \param[in] ref : the associated model instance.
             */
            Iterator (Model& ref)  : _ref(ref)   {}

            /** Set the data to be iterated.
             * \param[in] d : the data as information source for the iterator
             */
            void setData (tools::misc::Data& d);

        private:
            /** Reference on the underlying model; called for its 'build' method. */
            Model& _ref;
        };

        /************************************************************/
        /** \brief Iterator on neighbors kmers of a given kmer.
         *
         * This class will iterate the 8 neighbors kmers from a given source kmer.
         * It is similar to the Model::iterateNeighbors, except that here we don't have a container
         * holding all the successive kmers (ie. we have here only sequential access and
         * not direct access).
         *
         * To be used, such an iterator must be initialized with a source kmer,
         * which is done with the 'setSource' method.
         */
        class NeighborIterator : public tools::dp::impl::VectorIterator<Type>
        {
        public:

            /** Constructor.
             * \param[in] ref : the associated model instance. */
            NeighborIterator (Model& ref)  : _ref(ref) {  this->_items.resize(8);   this->_nb = 8; }

            /** Set the source kmer from which we want to iterate the neighbors.
             * \param[in] source : the source kmer. */
            void setSource (const Type& source);

        private:
            /** Reference on the underlying model; called for its 'build' method. */
            Model& _ref;
        };

    protected:

        /** Size of a kmer for this model. */
        size_t  _kmerSize;

        /** Mask for the kmer. Used for computing recursively kmers. */
        Type  _kmerMask;

        /** Shortcut for easing/speeding up the recursive revcomp computation. */
        Type _revcompTable[4];

        /** Mode for computing the kmers. */
        KmerMode _mode;

        /** Functor that returns the provided nucleotide.
         * \param[in] nt : the nucleotide in ASCII
         * \return the same value as input. */
        static int NTidentity(char nt)  { return nt; }

        /** Transform a nucleotide in ASCII form into an integer form as:
         *     - A=0
         *     - C=1
         *     - T=2
         *     - G=3
         * \param[in] nt : the nucleotide in ASCII
         * \return the translated nucleotide */
        static int NT2int(char nt)  {  return (nt>>1)&3;  }

        /** Generic function that computes a kmer from a sequence.
         * \param[in] seq : source data from which the kmer is computed.
         * \param[in] encoding : encoding of the source data
         * \param[in] mode : tells how to compute the kmer.
         * \return the computed kmer.
         */
        Type codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding, KmerMode mode)  const ;

        /** Compute the successor of a kmer in a recursive way.
         *  WARNING ! we can't compute the minimum of direct and revcomp with this function since we need to know
         *  both current direct and revcomp for computing the next one.
         * \param[in] seed : initial kmer from which we want to compute the successor
         * \param[in] nucl : the nucleotide to be appended to the current kmer
         * \param[in] encoding : encoding of the source data
         * \param[in] mode : tells how to compute the kmer.
         */
        Type codeSeedRight (const Type& seed, char nucl, tools::misc::Data::Encoding_e encoding, KmerMode mode)  const ;

        /** Build a kmer from a Data object (ie a sequence of nucleotides), starting at an index in the nucleotides sequence.
         * The result is a pair holding the built kmer and a boolean set to yes if the built kmer has to be understood in
         * the forward sense, false otherwise.
         * \param[in] data : the data from which we extract a kmer
         * \param[in] idx : start index in the data object (default to 0)
         * \param[in] mode : the mode for building kmers
         * \return a pair with the built kmer and a boolean set to yes if the kmer is understood in the forward strand
         */
        std::pair<Type, bool> getKmer (const tools::misc::Data& data, size_t idx, KmerMode mode)  const ;

        /** Build a vector of successive kmers from a given sequence of nucleotides provided as a Data object.
         * \param[in] data : the sequence of nucleotides.
         * \param[out] kmersBuffer : the successive kmers built from the data object.
         * \param[in] mode : the mode for building kmers
         * \return true if kmers have been extracted, false otherwise. */
        bool build (tools::misc::Data& data, std::vector<Type>& kmersBuffer, KmerMode mode)  const ;

        /** Iteration of the kmers from a data object through a functor (so lambda expressions can be used).
         * \param[in] data : the sequence of nucleotides.
         * \param[in] fct  : functor that handles one kmer
         * \return true if kmers have been extracted, false otherwise. */
        template<typename Functor> bool iterate (tools::misc::Data& data, Functor fct, KmerMode mode) const;

    };  // class Model


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

};  // class Kmer

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

/** We include template definitions. */
#include <gatb/kmer/impl/Model.tpp>

#endif /* _GATB_CORE_KMER_IMPL_MODEL_HPP_ */
