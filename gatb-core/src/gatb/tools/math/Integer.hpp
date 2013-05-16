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

#if INTEGER_KIND==1
    #include <gatb/tools/math/NativeInt64.hpp>
    namespace gatb  {  namespace core  { namespace tools {  namespace math  {
        typedef gatb::core::tools::math::NativeInt64 Integer;
    }}}};

#elif INTEGER_KIND==2
    #include <gatb/tools/math/NativeInt128.hpp>
    namespace gatb  {  namespace core  { namespace tools {  namespace math  {
        typedef gatb::core::tools::math::NativeInt128 Integer;
    }}}};

#else
    #include <gatb/tools/math/LargeInt.hpp>
    namespace gatb  {  namespace core  { namespace tools {  namespace math  {
        typedef gatb::core::tools::math::LargeInt<KMER_PRECISION> Integer;
    }}}};
#endif

/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_HPP_ */
