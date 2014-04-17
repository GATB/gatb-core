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

/** \file Integer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Entry point class for large integer usage
 *
 *  We define here (at compilation time) what will be the class used for large
 *  integers calculus.
 *
 *  According to the INTEGER_KIND compilation flag, we define the Integer class
 *  as an alias of one from several possible implementations.
 *
 *  Note that we have 2 possible native implementations (NativeInt64 and NativeInt128)
 *  that rely on native types uint64_t and __uint128_t.
 *
 *  For larger integer, a multi-precision LargeInt is used.
 *
 *  From the user point of view, [s]he has just to include this file and use the Integer
 *  class.
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_HPP_

/********************************************************************************/
#include <gatb/tools/math/LargeInt.hpp>

#include <boost/variant.hpp>

/********************************************************************************/
namespace gatb  {  namespace core  { namespace tools {  namespace math  {
/********************************************************************************/

template<typename T1, typename T2, typename T3, typename T4>
class IntegerTemplate
{
private:

    typedef boost::variant<T1,T2,T3,T4> Type;
    Type v;

          Type& operator *()       { return v; }
    const Type& operator *() const { return v; }

public:

    static char& getType()  {  static char instance = 0; return instance; }
    static void setType (char type)  {  getType() = type; }

    IntegerTemplate (int64_t n=0)
    {
        switch (getType())
        {
        case PREC_1: v = T1(n); break;
        case PREC_2: v = T2(n); break;
        case PREC_3: v = T3(n); break;
        case PREC_4: v = T4(n); break;
        default:  if (getType()<=PREC_4) { v = T4(n);  break; } else { throw "Integer not initialized"; }
        }
    }

    template<typename T>  explicit IntegerTemplate (const T& t) : v (t)  {}

    template<typename T>
    IntegerTemplate& operator=(const T& t)
    {
        v = t;
        return *this;
    }

    const char*  getName ()         { return boost::apply_visitor (Integer_name(),  *(*this)); }
    const size_t getSize ()         { return boost::apply_visitor (Integer_size(),  *(*this)); }
    hid_t hdf5 (bool& isCompound)   { return boost::apply_visitor (Integer_hdf5(isCompound),  *(*this)); }

    inline friend IntegerTemplate   operator+  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_plus(),  *a, *b);  }
    inline friend IntegerTemplate   operator-  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_minus(), *a, *b);  }
    inline friend IntegerTemplate   operator|  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_or(),    *a, *b);  }
    inline friend IntegerTemplate   operator^  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_xor(),   *a, *b);  }
    inline friend IntegerTemplate   operator&  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_and(),   *a, *b);  }
    inline friend IntegerTemplate   operator~  (const IntegerTemplate& a)                    {  Integer_compl v; return  boost::apply_visitor (v, *a);  }

    inline friend bool      operator== (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_equals(), *a, *b);  }
    inline friend bool      operator!= (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  ! (*a==*b);  }
    inline friend bool      operator<  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_less(),   *a, *b);  }
    inline friend bool      operator<= (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_lesseq(), *a, *b);  }

    inline friend IntegerTemplate   operator*  (const IntegerTemplate& a, const int&       c)  {  return  boost::apply_visitor (Integer_mult(c), *a);  }
    inline friend IntegerTemplate   operator/  (const IntegerTemplate& a, const u_int32_t& c)  {  return  boost::apply_visitor (Integer_div(c),  *a);  }
    inline friend u_int32_t operator%  (const IntegerTemplate& a, const u_int32_t& c)  {  return  boost::apply_visitor (Integer_mod(c),  *a);  }

    inline friend IntegerTemplate   operator>> (const IntegerTemplate& a, const int& c)  {  return  boost::apply_visitor (Integer_shiftLeft(c),   *a);  }
    inline friend IntegerTemplate   operator<< (const IntegerTemplate& a, const int& c)  {  return  boost::apply_visitor (Integer_shiftRight(c),  *a);  }

    void operator+= (const IntegerTemplate& a)  {  boost::apply_visitor (Integer_plusaffect(),  *(*this), *a);  }
    void operator^= (const IntegerTemplate& a)  {  boost::apply_visitor (Integer_xoraffect(),   *(*this), *a);  }

    u_int8_t  operator[]  (size_t idx) const   { return  boost::apply_visitor (Integer_value_at(idx), *(*this)); }

    friend IntegerTemplate revcomp (const IntegerTemplate& a,  size_t sizeKmer)  {  return  boost::apply_visitor (Integer_revomp(sizeKmer),  *a);  }

    friend u_int64_t hash1        (const IntegerTemplate& a,  u_int64_t seed)  {  return  boost::apply_visitor (Integer_hash1(seed),  *a);          }
    friend u_int64_t oahash       (const IntegerTemplate& a)                   {  return  boost::apply_visitor (Integer_oahash(), *a);              }
    friend u_int64_t simplehash16 (const IntegerTemplate& a,  int shift)       {  return  boost::apply_visitor (Integer_simplehash16(shift),  *a);  }

    std::string toString (size_t sizeKmer) const  {  return boost::apply_visitor (Integer_toString(sizeKmer), *(*this)); }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const IntegerTemplate& a)  {  s << *a;  return s;  }

    template<typename U>
    const U& get ()  const  {  return * boost::get<U>(&v);  }

private:

    struct Integer_name : public boost::static_visitor<const char*>    {
        template<typename T>  const char* operator() (const T& a) const { return a.getName();  }};

    struct Integer_size : public boost::static_visitor<const size_t>    {
        template<typename T>  const size_t operator() (const T& a) const  { return a.getSize();  }};

    struct Integer_plus : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a + b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();       }
    };

    struct Integer_minus : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a - b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_or : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a | b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_xor : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a ^ b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_and : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a & b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_less : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a < b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;  }
    };

    struct Integer_lesseq : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a <= b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;   }
    };

    struct Integer_equals : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a == b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;   }
    };

    struct Integer_plusaffect : public boost::static_visitor<>    {
        template<typename T>              void operator() ( T& a, const T& b) const  { a += b;  }
        template<typename T, typename U>  void operator() ( T& a, const U& b) const  {   }
    };

    struct Integer_xoraffect : public boost::static_visitor<>    {
        template<typename T>              void operator() ( T& a, const T& b) const  { a ^= b;  }
        template<typename T, typename U>  void operator() ( T& a, const U& b) const  {   }
    };

    struct Integer_compl : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>  IntegerTemplate operator() (const T& a)  { return IntegerTemplate(~a);  }};

    template<typename Result, typename Arg>
    struct Visitor : public boost::static_visitor<Result>
    {
        Visitor (Arg a=Arg()) : arg(a) {}
        Arg arg;
    };

    struct Integer_hdf5 : public Visitor<hid_t,bool&>   {
        Integer_hdf5 (bool& c) : Visitor<IntegerTemplate,bool&>(c) {}
        template<typename T>  hid_t operator() (const T& a)  { return a.hdf5 (this->arg);  }};

    struct Integer_mult : public Visitor<IntegerTemplate,const int>    {
        Integer_mult (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate(a*this->arg);  }};

    struct Integer_div : public Visitor<IntegerTemplate,const u_int32_t>    {
        Integer_div (const u_int32_t& c) : Visitor<IntegerTemplate,const u_int32_t>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate(a/this->arg);  }};

    struct Integer_mod : public Visitor<u_int32_t,const u_int32_t>    {
        Integer_mod (const u_int32_t& c) : Visitor<u_int32_t,const u_int32_t>(c) {}
        template<typename T>  u_int32_t operator() (const T& a) const  { return (a%this->arg);  }};

    struct Integer_shiftLeft : public Visitor<IntegerTemplate,const int>    {
        Integer_shiftLeft (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (a >> this->arg);  }};

    struct Integer_shiftRight : public Visitor<IntegerTemplate,const int>    {
        Integer_shiftRight (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (a << this->arg);  }};

    struct Integer_revomp : public Visitor<IntegerTemplate,size_t>    {
        Integer_revomp (const size_t& c) : Visitor<IntegerTemplate,size_t>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (revcomp(a,this->arg));  }};

    struct Integer_hash1 : public Visitor<u_int64_t,u_int64_t>    {
        Integer_hash1 (const u_int64_t& c) : Visitor<u_int64_t,u_int64_t>(c) {}
        template<typename T>  u_int64_t operator() (const T& a) const  { return (hash1(a,this->arg));  }};

    struct Integer_oahash : public boost::static_visitor<u_int64_t>    {
        template<typename T>  u_int64_t operator() (const T& a) const  { return (oahash(a));  }};

    struct Integer_simplehash16 : public Visitor<u_int64_t,int>    {
        Integer_simplehash16 (const int& c) : Visitor<u_int64_t,int>(c) {}
        template<typename T>  u_int64_t operator() (const T& a) const  { return (simplehash16(a,this->arg));  }};

    struct Integer_value_at : public Visitor<u_int8_t,size_t>   {
        Integer_value_at (size_t idx) : Visitor<u_int8_t,size_t>(idx) {}
        template<typename T>  u_int8_t operator() (const T& a) const { return a[this->arg];  }};

    struct Integer_toString : public Visitor<std::string,size_t>   {
        Integer_toString (size_t c) : Visitor<std::string,size_t>(c) {}
        template<typename T>  std::string operator() (const T& a) const  { return a.toString(this->arg);  }};
};

/********************************************************************************/

#define INTEGER_TYPES   LargeInt<PREC_1>,LargeInt<PREC_2>,LargeInt<PREC_3>,LargeInt<PREC_4>

typedef IntegerTemplate <INTEGER_TYPES> Integer;

/********************************************************************************/
}}}};
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_HPP_ */
