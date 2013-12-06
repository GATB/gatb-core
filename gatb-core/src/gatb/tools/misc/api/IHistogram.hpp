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

#include <gatb/system/api/ISmartPointer.hpp>
#include <hdf5/hdf5.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** Here is a command line for showing the histogram with gnuplot from the hdf5 file 'graph.h5'
    h5dump -y -d dsk/histogram graph.h5 | grep [0-9] | grep -v [A-Z].* | paste - - | gnuplot -p -e 'plot [][0:100] "-" with lines'

 For the sum of the distribution, you can use;
    h5dump -y -d dsk/histogram graph.h5 | grep [0-9] | grep -v [A-Z].* | paste - - | gawk 'BEGIN{s=0; i=0} { s=s+$2; i=i+1; print i,"  ", s}' | gnuplot -p -e 'plot [0:10][0:] "-" with lines'
*/

/** \brief TBD */
class IHistogram : virtual public system::ISmartPointer
{
public:

    /********************************************************************************/
    struct Entry
    {
        u_int16_t index;
        u_int64_t abundance;

        inline static hid_t hdf5 (bool& compound)
        {
            hid_t result = H5Tcreate (H5T_COMPOUND, sizeof(Entry));
            H5Tinsert (result, "index",      HOFFSET(Entry, index),     H5T_NATIVE_UINT16);
            H5Tinsert (result, "abundance",  HOFFSET(Entry, abundance), H5T_NATIVE_UINT64);
            compound = true;
            return result;
        }
    };

    /** */
    virtual ~IHistogram() {}

    /** */
    virtual size_t getLength() = 0;

    /** */
    virtual void inc (u_int16_t index) = 0;

    /** */
    virtual void save () = 0;

    /** */
    virtual u_int64_t& get (u_int16_t idx) = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_ */
