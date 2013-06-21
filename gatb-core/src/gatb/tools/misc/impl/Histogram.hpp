/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Histogram.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Histogram feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_
#define _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <string>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** TBD */
class Histogram
{
public:

    Histogram (size_t length, const std::string& uri) : _length(length), _uri(uri)  {  _histogram.resize (_length + 1); }

    void inc (size_t abundance);

    void save ();

private:

    size_t                 _length;
    std::string            _uri;
    std::vector<u_int64_t> _histogram;
};

/********************************************************************************/

/** */
class HistogramNull : public Histogram
{
public:
    HistogramNull () : Histogram (0,"") {}
    void inc ()  {}
    void save ()  {}
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_ */
