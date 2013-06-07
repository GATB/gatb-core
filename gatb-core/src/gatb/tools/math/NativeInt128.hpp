/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file NativeInt128.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native 128 bits integer type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_128_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_128_HPP_

/********************************************************************************/
#ifdef INT128_FOUND
/********************************************************************************/

#include <iostream>
#include <gatb/system/api/types.hpp>
#include <gatb/tools/math/NativeInt64.hpp>

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
class NativeInt128
{
public:

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    NativeInt128 (const __uint128_t& c=0)  {  value = c;  }

    static const char* getName ()  { return "NativeInt128"; }

    NativeInt128 operator+  (const NativeInt128& other)     const   {  return value + other.value;  }
    NativeInt128 operator-  (const NativeInt128& other)     const   {  return value - other.value;  }
    NativeInt128 operator*  (const int& coeff)              const   {  return value * coeff;        }
    NativeInt128 operator/  (const u_int32_t& divisor)      const   {  return value / divisor;      }
    u_int32_t    operator%  (const u_int32_t& divisor)      const   {  return value % divisor;      }
    NativeInt128 operator^  (const NativeInt128& other)     const   {  return value ^ other.value;  }
    NativeInt128 operator&  (const NativeInt128& other)     const   {  return value & other.value;  }
    NativeInt128 operator&  (const char& other)             const   {  return value & other;        }
    NativeInt128 operator~  ()                              const   {  return ~value;               }
    NativeInt128 operator<< (const int& coeff)              const   {  return value << coeff;       }
    NativeInt128 operator>> (const int& coeff)              const   {  return value >> coeff;       }
    bool         operator!= (const NativeInt128& c)         const   {  return value != c.value;     }
    bool         operator== (const NativeInt128& c)         const   {  return value == c.value;     }
    bool         operator<  (const NativeInt128& c)         const   {  return value < c.value;      }
    bool         operator<= (const NativeInt128& c)         const   {  return value <= c.value;     }

    NativeInt128& operator+=  (const NativeInt128& other)    {  value += other.value; return *this; }
    NativeInt128& operator^=  (const NativeInt128& other)    {  value ^= other.value; return *this; }

    /** Output stream overload. NOTE: for easier process, dump the value in hexadecimal.
     * \param[in] os : the output stream
     * \param[in] in : the integer value to be output.
     * \return the output stream.
     */
    friend std::ostream & operator<<(std::ostream & os, const NativeInt128 & in)
    {
        __uint128_t x = in.value;

        u_int64_t high_nucl = (u_int64_t) (x>>64);
        u_int64_t low_nucl  = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));

        if (high_nucl == 0) {   os << std::hex <<                     low_nucl << std::dec;  }
        else                {   os << std::hex << high_nucl << "." << low_nucl << std::dec;  }
        return os;
    }

    /********************************************************************************/
    
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=64).
     */
    inline void printASCII ( size_t sizeKmer = 64)
    {
        int i;
        u_int64_t temp = value;
        
        
        char seq[65];
        char bin2NT[4] = {'A','C','T','G'};
        
        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }
        seq[sizeKmer]='\0';
        
        std::cout << seq << std::endl;
    }
    
    
private:
    __uint128_t value;

    friend NativeInt128 revcomp (const NativeInt128& i,   size_t sizeKmer);
    friend u_int64_t    hash    (const NativeInt128& key, u_int64_t  seed);
};

/********************************************************************************/
inline NativeInt128 revcomp (const NativeInt128& in, size_t sizeKmer)
{
    //                  ---64bits--   ---64bits--
    // original kmer: [__high_nucl__|__low_nucl___]
    //
    // ex:            [         AC  | .......TG   ]
    //
    //revcomp:        [         CA  | .......GT   ]
    //                 \_low_nucl__/\high_nucl/

    const __uint128_t& x = in.value;

    u_int64_t high_nucl = (u_int64_t)(x>>64);
    int nb_high_nucl = sizeKmer>32?sizeKmer - 32:0;

    __uint128_t revcomp_high_nucl = NativeInt64::revcomp64 (high_nucl, nb_high_nucl);

    if (sizeKmer<=32) revcomp_high_nucl = 0; // srsly dunno why this is needed. gcc bug? u_int64_t x ---> (x>>64) != 0

    u_int64_t low_nucl = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));
    int nb_low_nucl = sizeKmer>32?32:sizeKmer;

    __uint128_t revcomp_low_nucl = NativeInt64::revcomp64 (low_nucl, nb_low_nucl);

    return (revcomp_low_nucl<<(2*nb_high_nucl)) + revcomp_high_nucl;
}

/********************************************************************************/
inline u_int64_t hash (const NativeInt128& item, u_int64_t seed)
{
    const __uint128_t& elem = item.value;

    return NativeInt64::hash64 ((u_int64_t)(elem>>64),seed) ^
           NativeInt64::hash64 ((u_int64_t)(elem&((((__uint128_t)1)<<64)-1)),seed);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

/********************************************************************************/
#endif //INT128_FOUND
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_128_HPP_ */
