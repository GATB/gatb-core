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
#include <gatb/tools/misc/api/Abundance.hpp>
#include <hdf5/hdf5.h>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

class NativeInt8 : private misc::ArrayData<u_int8_t, 1>
{
public:

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    NativeInt8 (const u_int8_t& c=0)  {  value[0] = c;  }

    static const char* getName ()  { return "NativeInt8"; }

    static size_t getSize ()  { return 8*sizeof(u_int8_t); }

    operator char ()  const {  return (char) value[0];  }

    NativeInt8 operator+  (const NativeInt8& other)   const   {  return value[0] + other.value[0];  }
    NativeInt8 operator-  (const NativeInt8& other)   const   {  return value[0] - other.value[0];  }
    NativeInt8 operator|  (const NativeInt8& other)   const   {  return value[0] | other.value[0];  }
    NativeInt8 operator^  (const NativeInt8& other)   const   {  return value[0] ^ other.value[0];  }
    NativeInt8 operator&  (const NativeInt8& other)   const   {  return value[0] & other.value[0];  }
    NativeInt8 operator&  (const char& other)          const   {  return value[0] & other;        }
    NativeInt8 operator~  ()                           const   {  return ~value[0];               }
    NativeInt8 operator<< (const int& coeff)           const   {  return value[0] << coeff;       }
    NativeInt8 operator>> (const int& coeff)           const   {  return value[0] >> coeff;       }
    bool        operator!= (const NativeInt8& c)       const   {  return value[0] != c.value[0];     }
    bool        operator== (const NativeInt8& c)       const   {  return value[0] == c.value[0];     }
    bool        operator<  (const NativeInt8& c)       const   {  return value[0] < c.value[0];      }
    bool        operator<= (const NativeInt8& c)       const   {  return value[0] <= c.value[0];     }
    NativeInt8& operator+=  (const NativeInt8& other)    {  value[0] += other.value[0]; return *this; }
    NativeInt8& operator^=  (const NativeInt8& other)    {  value[0] ^= other.value[0]; return *this; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const NativeInt8 & l)
    {
        s << std::hex << l.value[0] << std::dec;  return s;
    }

    /********************************************************************************/
    inline static hid_t hdf5 (bool& isCompound)
    {
        return H5Tcopy (H5T_NATIVE_UINT8);
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_8_HPP_ */
