/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
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
#include <hdf5.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for kmer management. */
namespace kmer      {
/********************************************************************************/

/** Define a kmer. */
template<typename Type, typename Number=u_int16_t> struct Kmer : Abundance<Type,Number>
{
    Kmer(const Type& val, const Number& abund) : Abundance<Type,Number>(val, abund) {}

    Kmer() : Abundance<Type,Number>(Type(), 0) {}

    Kmer(const Type& val) : Abundance<Type,Number>(val, 0) {}

    bool operator< (const Kmer& other) const {  return this->value < other.value; }


    /********************************************************************************/
    inline static hid_t hdf5 (bool& isCompound)
    {
        hid_t abundanceType = H5Tcopy (H5T_NATIVE_UINT16);

        hid_t result = H5Tcreate (H5T_COMPOUND, sizeof(Kmer));
        H5Tinsert (result, "value",      HOFFSET(Kmer, value),     Type::hdf5(isCompound));
        H5Tinsert (result, "abundance",  HOFFSET(Kmer, abundance), abundanceType);

        isCompound = true;

        return result;
    }
};

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

enum Strand
{
    STRAND_FORWARD = (1<<0),
    STRAND_REVCOMP = (1<<1),
    STRAND_ALL     = STRAND_FORWARD + STRAND_REVCOMP
};

inline Strand StrandReverse (const Strand& s)  {  return (s==STRAND_FORWARD ? STRAND_REVCOMP : STRAND_FORWARD);  }


enum Nucleotide
{
    NUCL_A   = 0,
    NUCL_C   = 1,
    NUCL_T   = 2,
    NUCL_G   = 3,
    NUCL_UNKNOWN = 4
};

inline char ascii (Nucleotide nt)
{
    static char table[] = {'A', 'C', 'T', 'G', 'N' };
    return table[(int)nt];
}

inline Nucleotide reverse (Nucleotide nt)
{
    static Nucleotide table[] = {NUCL_T, NUCL_G, NUCL_A, NUCL_C, NUCL_UNKNOWN};
    return table[(int)nt];
}

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
     * \param[in] encoding : encoding mode of the sequence
     * \return the kmer for the given nucleotides. */
    virtual kmer_type codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding) = 0;

    /** Compute the next right kmer given a current kmer and a nucleotide.
     * \param[in] val_seed : the current kmer as a starting point
     * \param[in] nucleotide : the next nucleotide
     * \param[in] encoding : encoding mode of the sequence
     * \param[in] mode : return mode of the kmer
     * \return the kmer on the right of the given kmer. */
    virtual kmer_type codeSeedRight (const kmer_type& val_seed, char nucleotide, tools::misc::Data::Encoding_e encoding, KmerMode mode = KMER_MINIMUM) = 0;

    /** Build a vector of kmers from some source data (as a Data instance).
     * According to the implementation of the IModel interface, the kmers can be different (direct, revcomp, minimum)
     * \param[in] d : source data from which the kmers are built
     * \param[in] kmersBuffer : vector where to put the built kmers into.
     * \param[in] mode : mode of the kmer enumeration
     */
    virtual bool build (tools::misc::Data& d, std::vector<kmer_type>& kmersBuffer, KmerMode mode = KMER_MINIMUM) = 0;

    /** Destructor. */
    virtual ~IModel () {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMODEL_HPP_ */
