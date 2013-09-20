/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file NativeInt8.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int8_t type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_8_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_8_HPP_

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
class NativeInt8
{
public:

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    NativeInt8 (const u_int8_t& c=0)  {  value = c;  }

    static const char* getName ()  { return "NativeInt8"; }

    static const size_t getSize ()  { return 8*sizeof(u_int8_t); }

    NativeInt8 operator+  (const NativeInt8& other)   const   {  return value + other.value;  }
    NativeInt8 operator-  (const NativeInt8& other)   const   {  return value - other.value;  }
    NativeInt8 operator|  (const NativeInt8& other)   const   {  return value | other.value;  }
    NativeInt8 operator^  (const NativeInt8& other)   const   {  return value ^ other.value;  }
    NativeInt8 operator&  (const NativeInt8& other)   const   {  return value & other.value;  }
    NativeInt8 operator&  (const char& other)          const   {  return value & other;        }
    NativeInt8 operator~  ()                           const   {  return ~value;               }
    NativeInt8 operator<< (const int& coeff)           const   {  return value << coeff;       }
    NativeInt8 operator>> (const int& coeff)           const   {  return value >> coeff;       }
    bool        operator!= (const NativeInt8& c)       const   {  return value != c.value;     }
    bool        operator== (const NativeInt8& c)       const   {  return value == c.value;     }
    bool        operator<  (const NativeInt8& c)       const   {  return value < c.value;      }
    bool        operator<= (const NativeInt8& c)       const   {  return value <= c.value;     }
    NativeInt8& operator+=  (const NativeInt8& other)    {  value += other.value; return *this; }
    NativeInt8& operator^=  (const NativeInt8& other)    {  value ^= other.value; return *this; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const NativeInt8 & l)
    {
        s << std::hex << l.value << std::dec;  return s;
    }

    /********************************************************************************/
    inline static hid_t hdf5 (bool& isCompound)
    {
        return H5Tcopy (H5T_NATIVE_UINT8);
    }
    
private:
    u_int8_t value;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_8_HPP_ */
