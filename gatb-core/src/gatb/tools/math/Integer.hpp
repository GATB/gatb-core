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

/** \file Integer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Smart Pointer Design Pattern definition
 *
 * arbitrary-precision integer library
 * very limited: only does what minia needs (but not what minia deserves)
 * This file holds interfaces related to the Design Pattern Observer.
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_HPP_

/********************************************************************************/

#include <stdint.h>
#include <algorithm>
#include <iostream>

#ifndef ASSERTS
#define NDEBUG // disable asserts; those asserts make sure that with PRECISION == [1 or 2], all is correct
#endif

// some 64-bit assert macros
#if defined(_LP64) && defined(_largeint)
#define assert128(x) assert(precision != 2 || (x));
#else
#define assert128(x) ;
#endif

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

//#define DYNAMIC

/** \brief Large integer class
 */
template<int nbbits>  class Integer
{
public:

    static const int precision = (nbbits+63)/64;

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    Integer(const uint64_t& c)
    {
#ifdef DYNAMIC
        array = new uint64_t [precision];
#endif
        array[0] = c;   for (int i = 1; i < precision; i++)  { array[i] = 0; }
    }

    /** Default constructor. */
    Integer()
    {
#ifdef DYNAMIC
        array = new uint64_t [precision];
#endif
    }

#ifdef DYNAMIC
    /** Copy constructor. */
    Integer (const Integer& c)
    {
        array = new uint64_t [precision];
        memcpy (array, c.array, sizeof(u_int64_t)*precision);
    }
#endif

#ifdef DYNAMIC
    /** Destructor */
    ~Integer ()  { delete[] array; }
#endif

    /********************************************************************************/
    Integer operator+(const Integer& other) const
    {
        Integer result;
        int carry = 0;
        for (int i = 0 ; i < precision ; i++)
        {
            result.array[i] = array[i] + other.array[i] + carry;
            carry = (result.array[i] < array[i]) ? 1 : 0;
        }

        assert(precision != 1 || (result == other.array[0] + array[0]));
        assert128(result.toInt128() == other.toInt128() + toInt128());
        return result;
    }

    /********************************************************************************/
    Integer operator-(const Integer& other) const
    {
        Integer result;
        int carry = 0;
        for (int i = 0 ; i < precision ; i++)
        {
            result.array[i] = array[i] - other.array[i] - carry;
            carry = (result.array[i] > array[i]) ? 1 : 0;
        }

        assert(precision != 1 || (result == array[0] - other.array[0]));
        assert128(result.toInt128() == toInt128() - other.toInt128());
        return result;
    }

    /********************************************************************************/
    Integer operator*(const int& coeff) const
    {
        Integer result (*this);
        // minia doesn't have that many multiplications cases

        if (coeff == 2 || coeff == 4)
        {
            result = result << (coeff / 2);
        }
        else
        {
            if (coeff == 21)
            {
                result = (result << 4) + (result << 2) + result;
            }
            else
            {
                printf("unsupported Integer multiplication: %d\n",coeff);
                exit(1);
            }
        }

        assert(precision != 1 || (result == array[0] * coeff));
        assert128(result.toInt128() == toInt128() * coeff);
        return result;
    }

    /********************************************************************************/
    Integer operator/(const uint32_t& divisor) const
    {
        Integer result;
        std::fill( result.array, result.array + precision, 0 );

        // inspired by Divide32() from http://subversion.assembla.com/svn/pxcode/RakNet/Source/BigInt.cpp

        uint64_t r = 0;
        uint32_t mask32bits = ~0;
        for (int i = precision-1; i >= 0; --i)
        {
            for (int j = 1; j >= 0; --j) // [j=1: high-32 bits, j=0: low-32 bits] of array[i]
            {
                uint64_t n = (r << 32) | ((array[i] >> (32*j)) & mask32bits );
                result.array[i] = result.array[i] | (((n / divisor) & mask32bits) << (32*j));
                r = n % divisor;
            }
        }
        assert(precision != 1 || (result == array[0] / divisor));
        assert128(result.toInt128() == toInt128() / divisor);
        return result;
    }

    /********************************************************************************/
    uint32_t operator%(const uint32_t& divisor) const
    {
        uint64_t r = 0;
        uint32_t mask32bits = ~0;
        for (int i = precision-1; i >= 0; --i)
        {
            for (int j = 1; j >= 0; --j) // [j=1: high-32 bits, j=0: low-32 bits] of array[i]
            {
                uint64_t n = (r << 32) | ((array[i] >> (32*j)) & mask32bits );
                r = n % divisor;
            }
        }

        assert(precision != 1 || (r == array[0] % divisor));
        assert128(r == toInt128() % divisor);
        return (uint32_t)r;
    }

    /********************************************************************************/
    Integer operator^(const Integer& other) const
    {
        Integer result;
        for (int i=0 ; i < precision ; i++)
            result.array[i] = array[i] ^ other.array[i];

        assert(precision != 1 || (result == (array[0] ^ other.array[0])));
        assert128(result.toInt128() == (toInt128() ^ other.toInt128()));
        return result;
    }

    /********************************************************************************/
    Integer operator&(const Integer& other) const
    {
        Integer result;
        for (int i=0 ; i < precision ; i++)
            result.array[i] = array[i] & other.array[i];

        assert(precision != 1 || (result == (array[0] & other.array[0])));
        assert128(result.toInt128() == (toInt128() & other.toInt128()));
        return result;
    }

    /********************************************************************************/
    Integer operator&(const char& other) const
    {
        Integer result;
        result.array[0] = array[0] & other;
        return result;
    }

    /********************************************************************************/
    Integer operator~() const
                    {
        Integer result;
        for (int i=0 ; i < precision ; i++)
            result.array[i] = ~array[i];

        assert(precision != 1 || (result == ~array[0]));
        assert128(result.toInt128() == ~toInt128());
        return result;
                    }

    /********************************************************************************/
    Integer operator<<(const int& coeff) const
    {
        Integer result (0);

        int large_shift = coeff / 64;
        int small_shift = coeff % 64;

        for (int i = large_shift ; i < precision-1; i++)
        {
            result.array[i] = result.array[i] | (array[i-large_shift] << small_shift);

            if (small_shift == 0) // gcc "bug".. uint64_t x; x>>64 == 1<<63, x<<64 == 1
            {
                result.array[i+1] = 0;
            }
            else
            {
                result.array[i+1] = array[i-large_shift] >> (64 - small_shift);
            }

        }
        result.array[precision-1] = result.array[precision-1] | (array[precision-1-large_shift] << small_shift);

        assert(precision != 1 || (result == (array[0] << coeff)));
        assert128(result.toInt128() == (toInt128() << coeff));
        return result;
    }

    /********************************************************************************/
    Integer operator>>(const int& coeff) const
    {
        Integer result (0);

        int large_shift = coeff / 64;
        int small_shift = coeff % 64;

        result.array[0] = (array[large_shift] >> small_shift);

        for (int i = 1 ; i < precision - large_shift ; i++)
        {
            result.array[i] = (array[i+large_shift] >> small_shift);
            if (small_shift == 0 && large_shift > 0) // gcc "bug".. uint64_t x; x>>64 == 1<<63, x<<64 == 1
            {
                result.array[i-1] =  result.array[i-1];
            }
            else
            {
                result.array[i-1] =  result.array[i-1] | (array[i+large_shift] << (64 - small_shift));
            }
        }

        assert(precision != 1 || ( small_shift == 0 || (result == array[0] >> coeff)));
        assert128(small_shift == 0 || (result.toInt128() == (toInt128() >> coeff)));
        return result;
    }

    /********************************************************************************/
    bool     operator!=(const Integer& c) const
                    {
        for (int i = 0 ; i < precision ; i++)
            if( array[i] != c.array[i] )
                return true;
        return false;
                    }

    /********************************************************************************/
    bool     operator==(const Integer& c) const
                    {
        for (int i = 0 ; i < precision ; i++)
            if( array[i] != c.array[i] )
                return false;
        return true;
                    }

    /********************************************************************************/
    bool     operator<(const Integer& c) const
    {
        for (int i = precision-1 ; i>=0 ; --i)
            if( array[i] != c.array[i] )
                return array[i] < c.array[i];

        return false;
    }

    /********************************************************************************/
    bool     operator<=(const Integer& c) const
                    {
        return operator==(c) || operator<(c);
                    }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const Integer & l)
    {
        s << l.toInt();  return s;
    }

    /********************************************************************************/
    uint64_t toInt() const
    {
        return array[0];
    }

    /********************************************************************************/
#ifdef _LP64
    __uint128_t toInt128() const
    {
        return ((__uint128_t)array[0]) + (((__uint128_t)array[1]) << ((__uint128_t)64));
    }
#endif

private:
#ifdef DYNAMIC
    uint64_t* array;
#else
    uint64_t array[precision];
#endif

    // c++ fun fact:
    // "const" will ban the function from being anything which can attempt to alter any member variables in the object.
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_HPP_ */
