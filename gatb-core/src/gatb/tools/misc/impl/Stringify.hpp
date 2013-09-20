/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Stringify.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool framework
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_STRINGIFY_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_STRINGIFY_HPP_

/********************************************************************************/

#include <string>
#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

class Stringify
{
public:

    static std::string format (const char* fmt, ...)
    {
        char* buffer = 0;

        va_list args;
        va_start (args, fmt);
        vasprintf (&buffer, fmt, args);
        va_end (args);

        if (buffer != NULL)
        {
            std::string result (buffer);
            free (buffer);
            return result;
        }
        else
        {
            return std::string ("");
        }
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_STRINGIFY_HPP_ */
