/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file NativeInt64.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_64_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_64_HPP_

/********************************************************************************/

#include <iostream>
#include <gatb/system/api/types.hpp>

extern const unsigned char revcomp_4NT[];
extern const unsigned char comp_NT    [];

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

/** \brief Large integer class
 */
class NativeInt64
{
public:

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    NativeInt64 (const u_int64_t& c=0)  {  value = c;  }

    static const char* getName ()  { return "NativeInt64"; }

    static const size_t getSize ()  { return 8*sizeof(value); }

    NativeInt64 operator+  (const NativeInt64& other)   const   {  return value + other.value;  }
    NativeInt64 operator-  (const NativeInt64& other)   const   {  return value - other.value;  }
    NativeInt64 operator*  (const int& coeff)           const   {  return value * coeff;        }
    NativeInt64 operator/  (const u_int32_t& divisor)   const   {  return value / divisor;      }
    u_int32_t   operator%  (const u_int32_t& divisor)   const   {  return value % divisor;      }
    NativeInt64 operator^  (const NativeInt64& other)   const   {  return value ^ other.value;  }
    NativeInt64 operator&  (const NativeInt64& other)   const   {  return value & other.value;  }
    NativeInt64 operator&  (const char& other)          const   {  return value & other;        }
    NativeInt64 operator~  ()                           const   {  return ~value;               }
    NativeInt64 operator<< (const int& coeff)           const   {  return value << coeff;       }
    NativeInt64 operator>> (const int& coeff)           const   {  return value >> coeff;       }
    bool        operator!= (const NativeInt64& c)       const   {  return value != c.value;     }
    bool        operator== (const NativeInt64& c)       const   {  return value == c.value;     }
    bool        operator<  (const NativeInt64& c)       const   {  return value < c.value;      }
    bool        operator<= (const NativeInt64& c)       const   {  return value <= c.value;     }

    NativeInt64& operator+=  (const NativeInt64& other)    {  value += other.value; return *this; }
    NativeInt64& operator^=  (const NativeInt64& other)    {  value ^= other.value; return *this; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const NativeInt64 & l)
    {
        s << std::hex << l.value << std::dec;  return s;
    }
    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=32).
     */
    inline void printASCII ( size_t sizeKmer = 32)
    {
        int i;
        u_int64_t temp = value;

        
        char seq[33];
        char bin2NT[4] = {'A','C','T','G'};
        
        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }
            seq[sizeKmer]='\0';
        
        std::cout << seq << std::endl;
    }
    
    /********************************************************************************/
    inline static u_int64_t revcomp64 (const u_int64_t& x, size_t sizeKmer)
    {
        u_int64_t res = x;

        unsigned char* kmerrev  = (unsigned char *) (&(res));
        unsigned char* kmer     = (unsigned char *) (&(x));

        for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        return (res >> (2*( 32 - sizeKmer))) ;
    }

    /********************************************************************************/
    inline static u_int64_t hash64 (u_int64_t key, u_int64_t seed)
    {
        u_int64_t hash = seed;
        hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
        hash = (~hash) + (hash << 21); // hash = (hash << 21) - hash - 1;
        hash = hash ^ (hash >> 24);
        hash = (hash + (hash << 3)) + (hash << 8); // hash * 265
        hash = hash ^ (hash >> 14);
        hash = (hash + (hash << 2)) + (hash << 4); // hash * 21
        hash = hash ^ (hash >> 28);
        hash = hash + (hash << 31);
        return hash;
    }

private:
    u_int64_t value;

    friend NativeInt64 revcomp (const NativeInt64& i,   size_t sizeKmer);
    friend u_int64_t    hash    (const NativeInt64& key, u_int64_t  seed);
};

/********************************************************************************/
inline NativeInt64 revcomp (const NativeInt64& x, size_t sizeKmer)
{
    return NativeInt64::revcomp64 (x.value, sizeKmer);
}

/********************************************************************************/
inline u_int64_t hash (const NativeInt64& key, u_int64_t seed)
{
    return NativeInt64::hash64 (key.value, seed);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_64_HPP_ */
