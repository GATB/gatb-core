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

/** \file LargeInt.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Class that manages large integers
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
#include <hdf5/hdf5.h>

#include <gatb/system/api/Exception.hpp>
#include <gatb/system/api/config.hpp>
#include <gatb/tools/math/NativeInt64.hpp>

#ifndef ASSERTS
#define NDEBUG // disable asserts; those asserts make sure that with PRECISION == [1 or 2], all is correct
#endif

#include <assert.h>

// some 64-bit assert macros
#if defined(_LP64)
#define assert128(x) assert(precision != 2 || (x));
#else
#define assert128(x) ;
#endif

extern const unsigned char revcomp_4NT[];
extern const unsigned char comp_NT    [];
#include <gatb/system/api/IMemory.hpp>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

/** \brief Large integer class
 *
 * The LargeInt class provides methods for integer calculus. It has a template parameter
 * 'precision' giving the number of bits used the integer representation. For instance:
 *  - LargeInt<1>  : representation of integers up to 2^64
 *  - LargeInt<2>  : representation of integers up to 2^128
 *  - etc
 *
 *  This template class has a specialization for precision=1. In this case, native 64 bits
 *  integers are used.
 *
 *  This template class may have a specialization for precision=2. If the used operating
 *  system allows it, native 128 bits integers are used.
 *
 *  In the other cases, the LargeInt provides a generic integer calculus class. Note that
 *  such an implementation could be optimized in several ways, including direct assembly
 *  code for maximum speed.
 *
 *  The LargeInt class is hugely used throughout the GATB project since it encodes kmers values.
 *
 *  The LargeInt class is mainly used with the IntegerTemplate class, where 4 specializations
 *  of LargeInt are used as template types of IntegerTemplate.
 *
 *  \see IntegerTemplate
 */
template<int precision>  class LargeInt : public misc::ArrayData<u_int64_t, precision>
{
public:

    /** Get the name of the class used by the variant (ie. one of the Ti template class parameters)
     * \return the class name.
     */
    static const char* getName ()
    {
        static char buffer[256];
        static bool first = true;
        if (first)  {  first = false;  snprintf (buffer, sizeof(buffer), "LargeInt<%d>", precision);  }
        return buffer;
    }

    /** Get the 64 less significant bits of the LargeInt object as a native integer type.
     * \return (part of) the LargeInt object as a native integer type.
     */
    u_int64_t getVal() const  { return this->value[0]; }

    /** Get the size of an instance of the class
     * \return the size of an object (in bits).
     */
    static const size_t getSize ()  { return 8*sizeof(u_int64_t)*precision; }

    /********************************************************************************/
    /** Constructor.
     * \param[in] val : initial value of the large integer. */
    LargeInt(const u_int64_t& val = 0)
    {
        this->value[0] = val;   for (int i = 1; i < precision; i++)  {  this->value[i] = 0;  }
    }

    /********************************************************************************/
    /** Returns lower 64 bits */
    u_int64_t toInt () const  {  throw system::Exception ("LargeInt<%d> no support of toInt", precision);  }

    /********************************************************************************/
    /** Operator +
     * \param[in] other : operand
     * \return sum of object and the operand.
     */
    LargeInt operator+ (const LargeInt& other) const
    {
        LargeInt result;
        int carry = 0;
        for (int i = 0 ; i < precision ; i++)
        {
            result.value[i] = this->value[i] + other.value[i] + carry;
            carry = (result.value[i] < this->value[i]) ? 1 : 0;
        }

        assert    (precision != 1 || (result == other.value[0] + this->value[0]));
        assert128 (result.toInt128() == other.toInt128() + toInt128());
        return result;
    }

    /********************************************************************************/
    /** Operator -
     * \param[in] other : operand
     * \return subtraction of object and the operand.
     */
    LargeInt operator- (const LargeInt& other) const
    {
        LargeInt result;
        int carry = 0;
        for (int i = 0 ; i < precision ; i++)
        {
            result.value[i] = this->value[i] - other.value[i] - carry;
            carry = (result.value[i] > this->value[i]) ? 1 : 0;
        }

        assert(precision != 1 || (result == this->value[0] - other.value[0]));
        assert128(result.toInt128() == toInt128() - other.toInt128());
        return result;
    }

    /********************************************************************************/
    /** Operator *
     * \param[in] coeff : operand
     * \return multiplication of the object by the coeff.
     */
    LargeInt operator*(const int& coeff) const
    {
        LargeInt result (*this);
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
                printf("unsupported LargeInt multiplication: %d\n",coeff);
                exit(1);
            }
        }

        assert(precision != 1 || (result == this->value[0] * coeff));
        assert128(result.toInt128() == toInt128() * coeff);
        return result;
    }

    /********************************************************************************/
    /** Operator /
     * \param[in] divisor : operand
     * \return division of the object by the divisor.
     */
    LargeInt operator/(const uint32_t& divisor) const
    {
        LargeInt result;
        std::fill( result.value, result.value + precision, 0 );

        // inspired by Divide32() from http://subversion.assembla.com/svn/pxcode/RakNet/Source/BigInt.cpp

        u_int64_t r = 0;
        uint32_t mask32bits = ~0;
        for (int i = precision-1; i >= 0; --i)
        {
            for (int j = 1; j >= 0; --j) // [j=1: high-32 bits, j=0: low-32 bits] of array[i]
            {
                u_int64_t n = (r << 32) | ((this->value[i] >> (32*j)) & mask32bits );
                result.value[i] = result.value[i] | (((n / divisor) & mask32bits) << (32*j));
                r = n % divisor;
            }
        }
        assert(precision != 1 || (result == this->value[0] / divisor));
        assert128(result.toInt128() == toInt128() / divisor);
        return result;
    }

    /********************************************************************************/
    /** Operator %
     * \param[in] divisor: operand
     * \return modulo of the object by the operand.
     */
    uint32_t operator%(const uint32_t& divisor) const
    {
        u_int64_t r = 0;
        uint32_t mask32bits = ~0;
        for (int i = precision-1; i >= 0; --i)
        {
            for (int j = 1; j >= 0; --j) // [j=1: high-32 bits, j=0: low-32 bits] of array[i]
            {
                u_int64_t n = (r << 32) | ((this->value[i] >> (32*j)) & mask32bits );
                r = n % divisor;
            }
        }

        assert(precision != 1 || (r == this->value[0] % divisor));
        assert128(r == toInt128() % divisor);
        return (uint32_t)r;
    }

    /********************************************************************************/
    /** Operator ^
     * \param[in] other : operand
     * \return operator^ of the object by the operand.
     */
    LargeInt operator^(const LargeInt& other) const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = this->value[i] ^ other.value[i];

        assert(precision != 1 || (result == (this->value[0] ^ other.value[0])));
        assert128(result.toInt128() == (toInt128() ^ other.toInt128()));
        return result;
    }

    /********************************************************************************/
    /** Operator |
     * \param[in] other : operand
     * \return operator| of the object by the operand.
     */
    LargeInt operator|(const LargeInt& other) const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = this->value[i] | other.value[i];
        
        assert(precision != 1 || (result == (this->value[0] | other.value[0])));
        assert128(result.toInt128() == (toInt128() | other.toInt128()));
        return result;
    }
    
    /********************************************************************************/
    /** Operator &
     * \param[in] other : operand
     * \return operator& of the object by the operand.
     */
    LargeInt operator&(const LargeInt& other) const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = this->value[i] & other.value[i];

        assert(precision != 1 || (result == (this->value[0] & other.value[0])));
        assert128(result.toInt128() == (toInt128() & other.toInt128()));
        return result;
    }

    /********************************************************************************/
    /** Operator &
     * \param[in] other : operand
     * \return operator& of the object by the operand.
     */
    LargeInt operator&(const char& other) const
    {
        LargeInt result;
        result.value[0] = this->value[0] & other;
        return result;
    }

    /********************************************************************************/
    /** Operator ~
     * \return negation of the object
     */
    LargeInt operator~() const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = ~this->value[i];

        assert(precision != 1 || (result == ~this->value[0]));
        assert128(result.toInt128() == ~toInt128());
        return result;
    }

    /********************************************************************************/
    /** Operator <<. Note: this method is likely to be hugely used when we want to get
     * neighbors of a given kmer encoded as a LargeInt object.
     * \param[in] coeff : operand
     * \return left shift of the object
     */
    LargeInt operator<<(const int& coeff) const
    {
        LargeInt result (0);

        int large_shift = coeff / 64;
        int small_shift = coeff % 64;

        for (int i = large_shift ; i < precision-1; i++)
        {
            result.value[i] = result.value[i] | (this->value[i-large_shift] << small_shift);

            if (small_shift == 0) // gcc "bug".. u_int64_t x; x>>64 == 1<<63, x<<64 == 1
            {
                result.value[i+1] = 0;
            }
            else
            {
                result.value[i+1] = this->value[i-large_shift] >> (64 - small_shift);
            }

        }
        result.value[precision-1] = result.value[precision-1] | (this->value[precision-1-large_shift] << small_shift);

        assert(precision != 1 || (result == (this->value[0] << coeff)));
        assert128(result.toInt128() == (toInt128() << coeff));
        return result;
    }

    /********************************************************************************/
    /** Operator >>. Note: this method is likely to be hugely used when we want to get
     * neighbors of a given kmer encoded as a LargeInt object.
     * \param[in] coeff : operand
     * \return right shift of the object
     */
    LargeInt operator>>(const int& coeff) const
    {
        LargeInt result (0);

        int large_shift = coeff / 64;
        int small_shift = coeff % 64;

        result.value[0] = (this->value[large_shift] >> small_shift);

        for (int i = 1 ; i < precision - large_shift ; i++)
        {
            result.value[i] = (this->value[i+large_shift] >> small_shift);
            if (small_shift == 0) // gcc "bug".. u_int64_t x; x>>64 == 1<<63, x<<64 == 1
            {
                result.value[i-1] =  result.value[i-1];
            }
            else
            {
                result.value[i-1] =  result.value[i-1] | (this->value[i+large_shift] << (64 - small_shift));
            }
        }

        assert(precision != 1 || ( small_shift == 0 || (result == this->value[0] >> coeff)));
        assert128(small_shift == 0 || (result.toInt128() == (toInt128() >> coeff)));
        return result;
    }

    /********************************************************************************/
    /** Operator !=
     * \param[in] c: operand
     * \return inequality
     */
    bool operator!=(const LargeInt& c) const
    {
        for (int i = 0 ; i < precision ; i++)
            if( this->value[i] != c.value[i] )
                return true;
        return false;
    }

    /********************************************************************************/
    /** Operator ==
     * \param[in] c: operand
     * \return equality
     */
    bool operator==(const LargeInt& c) const
    {
        for (int i = 0 ; i < precision ; i++)
            if( this->value[i] != c.value[i] )
                return false;
        return true;
    }

    /********************************************************************************/
    /** Operator <
     * \param[in] c: operand
     */
    bool operator<(const LargeInt& c) const
    {
        for (int i = precision-1 ; i>=0 ; --i)
            if( this->value[i] != c.value[i] )
                return this->value[i] < c.value[i];

        return false;
    }

    /********************************************************************************/
    /** Operator <=
     * \param[in] c: operand
     */
    bool operator<=(const LargeInt& c) const
    {
        return operator==(c) || operator<(c);
    }

    /********************************************************************************/
    /** Operator +=
     * \param[in] other : operand
     * \return addition and affectation
     */
    LargeInt& operator+=  (const LargeInt& other)
    {
        // NOT so easy to optimize because of the carry
        *this = *this + other;
        return *this;
    }

    /********************************************************************************/
    /** Operator ^=
     * \param[in] other : operand
     * \return xor and affectation
     */
    LargeInt& operator^=  (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  this->value[i] ^= other.value[i];  }
        return *this;
    }

    /********************************************************************************/
    /** Operator &=
     * \param[in] other : operand
     * \return and and affectation
     */
    LargeInt& operator&=  (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  this->value[i] &= other.value[i];  }
        return *this;
    }

    /********************************************************************************/
    /** Operator |=
     * \param[in] other : operand
     * \return or and affectation
     */
    LargeInt& operator|=  (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  this->value[i] |= other.value[i];  }
        return *this;
    }

    /********************************************************************************/
    /** Operator <<=
     * \param[in] coeff : operand
     * \return left shift and affectation
     */
    LargeInt& operator<<=  (const int& coeff)
    {
        *(this) = (*this) << coeff;  return *this;
    }

    /********************************************************************************/
    /** Operator >>=
     * \param[in] coeff : operand
     * \return right shift and affectation
     */
    LargeInt& operator>>=  (const int& coeff)
    {
        *(this) = (*this) >> coeff;  return *this;
    }

    /********************************************************************************/
    /** Synchronized atomic operation
     * \param[in] other : operand
     */
    LargeInt& sync_fetch_and_or (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  __sync_fetch_and_or (this->value + i, other.value[i]); }
        return *this;
    }

    /********************************************************************************/
    /** Synchronized atomic operation
     * \param[in] other : operand
     */
    LargeInt& sync_fetch_and_and (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  __sync_fetch_and_and (this->value + i, other.value[i]); }
        return *this;
    }

    /********************************************************************************/
    /** Output stream operator for the IntegerTemplate class
     * \param[in] s : the output stream to be used.
     * \param[in] l : the object to output
     * \return the modified output stream.
     */
    friend std::ostream & operator<<(std::ostream & s, const LargeInt<precision> & l)
    {
        int i=0;

        /** We want to display the number in hexa (easier to do...) */
        s << std::hex;

        /** We skip leading 0. */
        for (i=precision-1; i>=0 && l.value[i]==0; i--)  {}

        /** We dump the different parts of the large integer. */
        for (  ; i>=0 ; i--)  { s << l.value[i];   if (i>=1) { s << ".";  }  }

        /** We go back to decimal format. */
        s << std::dec;

        /** We return the output stream. */
        return s;
    }

    /********************************************************************************/
    /** Get the HDF5 type for the the class used by the variant  (ie. one of the Ti template class parameters)
     * \param[in] isCompound : tells whether the type is composed or not
     * \return a HDF5 identifier for the type.
     */
    inline static hid_t hdf5 (bool& isCompound)
    {
        hid_t result = H5Tcopy (H5T_NATIVE_INT);
        H5Tset_precision (result, 64*precision);
        return result;
    }

    /********************************************************************************/
    /** Computes a kmer value as polynom. We may have conversion from the data buffer to
     * a nucleotide code. This is done through the provided functor.
     * \param[in] data : kmer given as a buffer of nucleotides
     * \param[in] size : size of the kmer
     * \param[in] fct : convert the ith entry in the buffer into a nucleotide code  (A=0, C=1, T=2 and G=3)
     */
    template<typename Map>
    static LargeInt polynom (const char* data, size_t size, Map fct)
    {
        LargeInt res (0);
        for (size_t i=0; i<size; ++i)  {  res = res * 4 + fct(data[i]);  }
        return res;
    }

    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[in] sizeKmer : kmer size (def=32).
     * \return the ASCII string
     */
    std::string toString (size_t sizeKmer) const
    {
        char seq[sizeKmer+1];
        char bin2NT[4] = {'A','C','T','G'};

        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }
        seq[sizeKmer]='\0';
        return seq;
    }

    /********************************************************************************/
    /** Operator[] access the ith nucleotide in the given integer. For instance a[4] get the 5th nucleotide of
     * a kmer encoded as an Integer object.
     * \param[in] idx : index of the nucleotide to be retrieved
     * \return the nucleotide value as follow: A=0, C=1, T=2 and G=3
     */
    u_int8_t  operator[]  (size_t idx) const    {  return (this->value[idx/32] >> (2*idx)) & 3; }

private:

    template<int T>  friend LargeInt<T> revcomp (const LargeInt<T>& i, size_t sizeKmer);
    template<int T>  friend u_int64_t   hash1    (const LargeInt<T>& key, u_int64_t  seed);
    template<int T>  friend u_int64_t   oahash  (const LargeInt<T>& key);
    template<int T>  friend u_int64_t   simplehash16    (const LargeInt<T>& key, int  shift);

    // c++ fun fact:
    // "const" will ban the function from being anything which can attempt to alter any member variables in the object.
    // 2 months later: I don't understand this fun fact anymore. thanks c++.
};

/********************************************************************************/
template<int precision>  inline LargeInt<precision> revcomp (const LargeInt<precision>& x, size_t sizeKmer)
{
    const LargeInt<precision> res = x;

    unsigned char* kmerrev  = (unsigned char *) (&(res.value[0]));
    unsigned char* kmer     = (unsigned char *) (&(x.value[0]));

    for (size_t i=0; i<8*precision; ++i)
    {
        kmerrev[8*precision-1-i] = revcomp_4NT [kmer[i]];
    }

    return (res >> (2*( 32*precision - sizeKmer))  ) ;
}

/********************************************************************************/
template<int precision>  inline u_int64_t hash1 (const LargeInt<precision>& elem, u_int64_t seed=0)
{
    // hash = XOR_of_series[hash(i-th chunk iof 64 bits)]
    u_int64_t result = 0, chunk, mask = ~0;

    LargeInt<precision> intermediate = elem;
    for (size_t i=0;i<precision;i++)
    {
        chunk = (intermediate & mask).value[0];
        intermediate = intermediate >> 64;

        result ^= NativeInt64::hash64 (chunk,seed);
    }
    return result;
}

/********************************************************************************/
template<int precision>  u_int64_t oahash (const LargeInt<precision>& elem)
{
    // hash = XOR_of_series[hash(i-th chunk iof 64 bits)]
    u_int64_t result = 0, chunk, mask = ~0;

    LargeInt<precision> intermediate = elem;
    for (size_t i=0;i<precision;i++)
    {
        chunk = (intermediate & mask).value[0];
        intermediate = intermediate >> 64;
        result ^= NativeInt64::oahash64 (chunk);
    }
    return result;
}

	

/********************************************************************************/
template<int precision>  inline u_int64_t simplehash16 (const LargeInt<precision>& elem, int  shift)
{
    u_int64_t result = 0, chunk, mask = ~0;
    LargeInt<precision> intermediate = elem;
    
    chunk = (intermediate & mask).value[0];
    result ^= NativeInt64::simplehash16_64 (chunk,shift);

    return result;
}


/********************************************************************************/
/********************     SPECIALIZATION FOR precision=1     ********************/
/********************************************************************************/
#include <gatb/tools/math/LargeInt1.pri>

/********************************************************************************/
/********************     SPECIALIZATION FOR precision=2     ********************/
/********************************************************************************/
#include <gatb/tools/math/LargeInt2.pri>

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_ */
