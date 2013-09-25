/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file NativeInt16.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int16_t type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_16_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_16_HPP_

/********************************************************************************/

#include <iostream>
#include <gatb/system/api/types.hpp>
#include <hdf5.h>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

/** \brief Large integer class
 */
class NativeInt16 : private ArrayData<u_int16_t, 1>
{
public:

    typedef ArrayData<u_int16_t, 1> POD;

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    NativeInt16 (const u_int8_t& c=0)  {  value[0] = c;  }

    static const char* getName ()  { return "NativeInt16"; }

    static const size_t getSize ()  { return 8*sizeof(u_int16_t); }

    NativeInt16 operator+  (const NativeInt16& other)   const   {  return value[0] + other.value[0];  }
    NativeInt16 operator-  (const NativeInt16& other)   const   {  return value[0] - other.value[0];  }
    NativeInt16 operator|  (const NativeInt16& other)   const   {  return value[0] | other.value[0];  }
    NativeInt16 operator^  (const NativeInt16& other)   const   {  return value[0] ^ other.value[0];  }
    NativeInt16 operator&  (const NativeInt16& other)   const   {  return value[0] & other.value[0];  }
    NativeInt16 operator&  (const char& other)          const   {  return value[0] & other;        }
    NativeInt16 operator~  ()                           const   {  return ~value[0];               }
    NativeInt16 operator<< (const int& coeff)           const   {  return value[0] << coeff;       }
    NativeInt16 operator>> (const int& coeff)           const   {  return value[0] >> coeff;       }
    bool        operator!= (const NativeInt16& c)       const   {  return value[0] != c.value[0];     }
    bool        operator== (const NativeInt16& c)       const   {  return value[0] == c.value[0];     }
    bool        operator<  (const NativeInt16& c)       const   {  return value[0] < c.value[0];      }
    bool        operator<= (const NativeInt16& c)       const   {  return value[0] <= c.value[0];     }
    bool        operator>= (const NativeInt16& c)       const   {  return value[0] >= c.value[0];     }
    NativeInt16& operator+=  (const NativeInt16& other)    {  value[0] += other.value[0]; return *this; }
    NativeInt16& operator^=  (const NativeInt16& other)    {  value[0] ^= other.value[0]; return *this; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const NativeInt16 & l)
    {
        s << std::hex << l.value[0] << std::dec;  return s;
    }

    /********************************************************************************/
    inline static hid_t hdf5 (bool& isCompound)
    {
        return H5Tcopy (H5T_NATIVE_UINT16);
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_16_HPP_ */
