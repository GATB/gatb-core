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
namespace gatb      {
namespace core      {
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

enum Strand
{
    STRAND_FORWARD = (1<<0),
    STRAND_REVCOMP = (1<<1),
    STRAND_ALL     = STRAND_FORWARD + STRAND_REVCOMP
};

inline std::string toString (Strand s)
{
         if (s == STRAND_FORWARD)  { return std::string("FWD"); }
    else if (s == STRAND_REVCOMP)  { return std::string("REV"); }
    else { return std::string("???"); }
}

inline Strand StrandReverse (const Strand& s)  {  return (s==STRAND_FORWARD ? STRAND_REVCOMP : STRAND_FORWARD);  }

/********************************************************************************/

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
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMODEL_HPP_ */
