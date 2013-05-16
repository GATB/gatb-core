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

/** \file NativeInt128.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native 128 bits integer type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_128_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_128_HPP_

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
    NativeInt128 (const u_int128_t& c=0)  {  value = c;  }

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

    /** Output stream overload. NOTE: for easier process, dump the value in hexadecimal.
     * \param[in] os : the output stream
     * \param[in] in : the integer value to be output.
     * \return the output stream.
     */
    friend std::ostream & operator<<(std::ostream & os, const NativeInt128 & in)
    {
        u_int128_t x = in.value;

        u_int64_t high_nucl = (u_int64_t) (x>>64);
        u_int64_t low_nucl  = (u_int64_t)(x&((((u_int128_t)1)<<64)-1));

        if (high_nucl == 0) {   os << std::hex <<                     low_nucl << std::dec;  }
        else                {   os << std::hex << high_nucl << "." << low_nucl << std::dec;  }
        return os;
    }

private:
    u_int128_t value;

    friend NativeInt128 revcomp (NativeInt128& i, size_t sizeKmer);
};

/********************************************************************************/
inline NativeInt128 revcomp (NativeInt128& in, size_t sizeKmer)
{
    //                  ---64bits--   ---64bits--
    // original kmer: [__high_nucl__|__low_nucl___]
    //
    // ex:            [         AC  | .......TG   ]
    //
    //revcomp:        [         CA  | .......GT   ]
    //                 \_low_nucl__/\high_nucl/

    u_int128_t& x = in.value;

    u_int64_t high_nucl = (u_int64_t)(x>>64);
    int nb_high_nucl = sizeKmer>32?sizeKmer - 32:0;

    u_int128_t revcomp_high_nucl = NativeInt64::revcomp64 (high_nucl, nb_high_nucl);

    if (sizeKmer<=32) revcomp_high_nucl = 0; // srsly dunno why this is needed. gcc bug? u_int64_t x ---> (x>>64) != 0

    u_int64_t low_nucl = (u_int64_t)(x&((((u_int128_t)1)<<64)-1));
    int nb_low_nucl = sizeKmer>32?32:sizeKmer;

    u_int128_t revcomp_low_nucl = NativeInt64::revcomp64 (low_nucl, nb_low_nucl);

    return (revcomp_low_nucl<<(2*nb_high_nucl)) + revcomp_high_nucl;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_128_HPP_ */
