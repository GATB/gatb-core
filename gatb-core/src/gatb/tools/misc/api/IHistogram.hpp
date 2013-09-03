/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file IHistogram.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface for histogram (something counting abundances).
 */

#ifndef _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_
#define _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_

#include <gatb/tools/designpattern/api/SmartPointer.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief TBD */
class IHistogram : public virtual dp::ISmartPointer
{
public:

    /** */
    struct Entry
    {
        u_int16_t index;
        u_int32_t abundance;
    };

    /** */
    virtual size_t getLength() = 0;

    /** */
    virtual void inc (u_int32_t abundance) = 0;

    /** */
    virtual void save () = 0;

    /** */
    virtual u_int32_t& get (u_int16_t idx) = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_ */
