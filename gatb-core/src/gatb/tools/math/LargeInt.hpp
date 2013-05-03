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

/** \file LargeInt.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Smart Pointer Design Pattern definition
 *
 * arbitrary-precision integer library
 * very limited: only does what minia needs (but not what minia deserves)
 * This file holds interfaces related to the Design Pattern Observer.
 */

#ifndef _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_
#define _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_

/********************************************************************************/

#include <stdint.h>
#include <algorithm>
#include <iostream>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

/** \brief Large integer class
 */
template<int precision>  class LargeInt
{
public:

    /** Constructor.
     * \param[in] val : initial value of the large integer. */
    LargeInt(const uint64_t& val);
    LargeInt();

    // overloading
    LargeInt operator+(const LargeInt &) const;
    LargeInt operator-(const LargeInt &) const;
    LargeInt operator*(const int &) const;
    LargeInt operator/(const uint32_t &) const;
    uint32_t operator%(const uint32_t &) const;
    LargeInt operator^(const LargeInt &) const;
    LargeInt operator&(const LargeInt &) const;
    LargeInt operator~() const;
    LargeInt operator<<(const int &) const;
    LargeInt operator>>(const int &) const;
    bool     operator!=(const LargeInt &) const;
    bool     operator==(const LargeInt &) const;
    bool     operator<(const LargeInt &) const;
    bool     operator<=(const LargeInt &) const;

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const LargeInt<precision> & l)
    {
        s << l.toInt();  return s;
    }

    // custom
    uint64_t toInt() const;
    #ifdef _LP64
    __uint128_t toInt128() const;
    #endif

private:
    uint64_t array[precision];

    // c++ fun fact:
    // "const" will ban the function from being anything which can attempt to alter any member variables in the object.
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_ */
