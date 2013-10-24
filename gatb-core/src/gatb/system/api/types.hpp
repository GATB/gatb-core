/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file types.hpp
 *  \brief types definition for PLAST.
 *  \date 01/03/2013
 *  \author edrezen
 *
 *   We define here some types used thoughout the code.
 *
 *   Important: we define typedefs such as int16_t or u_int64_t. It is a good idea to use such typedefs
 *   instead of direct 'unsigned long' or 'short' for instance, because the actual number of used bytes
 *   may depend on the operating system/architecture. Using u_int32_t for instance ensure that we get
 *   an unsigned integer on 4 bytes.
 *
 *   Note that we use the <sys/types.h> file on Linux and MacOs. Such file may not exist on Windows (on Mingw
 *   to be more precise), so we propose here a definition. This is not perfect and should be improved.
 */

/********************************************************************************/

#ifndef _GATB_CORE_SYSTEM_TYPES_HPP_
#define _GATB_CORE_SYSTEM_TYPES_HPP_

/********************************************************************************/

#include <sys/types.h>

/********************************************************************************/

template<typename Type, int precision>
struct ArrayData
{
    Type value[precision];
};

/********************************************************************************/

/** Define an abundance. */
template<typename Type, typename Number=u_int16_t> struct Abundance
{
    Abundance (const Type& val=0, const Number& abund=0) : value(val), abundance(abund) {}

    Abundance& operator=(const Abundance& a)
    {
        if (&a != this)  {  value = a.value;  abundance=a.abundance;  }
        return *this;
    }

    const Number& getAbundance() const { return abundance; }
    const Type&   getValue()     const { return value;     }

    Type    value;
    Number  abundance;
};

/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_TYPES_HPP_ */
