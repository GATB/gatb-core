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

/** \file NativeInt64.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native uint64_t type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_64_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_64_HPP_

/********************************************************************************/

#include <iostream>
#include <stdint.h>

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
    NativeInt64 (const uint64_t& c=0)  {  value = c;  }

    static const char* getName ()  { return "NativeInt64"; }

    NativeInt64 operator+  (const NativeInt64& other)   const   {  return value + other.value;  }
    NativeInt64 operator-  (const NativeInt64& other)   const   {  return value - other.value;  }
    NativeInt64 operator*  (const int& coeff)           const   {  return value * coeff;        }
    NativeInt64 operator/  (const uint32_t& divisor)    const   {  return value / divisor;      }
    uint32_t    operator%  (const uint32_t& divisor)    const   {  return value % divisor;      }
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

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const NativeInt64 & l)
    {
        s << std::hex << l.value << std::dec;  return s;
    }

    /********************************************************************************/
    inline static uint64_t revcomp64 (uint64_t& x, size_t sizeKmer)
    {
        uint64_t res = x;

        unsigned char* kmerrev  = (unsigned char *) (&(res));
        unsigned char* kmer     = (unsigned char *) (&(x));

        for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        return (res >> (2*( 32 - sizeKmer))) ;
    }

private:
    uint64_t value;

    friend NativeInt64 revcomp (NativeInt64& i, size_t sizeKmer);
};

/********************************************************************************/
inline NativeInt64 revcomp (NativeInt64& x, size_t sizeKmer)
{
    return NativeInt64::revcomp64 (x.value, sizeKmer);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_64_HPP_ */
