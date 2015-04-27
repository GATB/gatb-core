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

/** \brief String helper for using printf-like formats
 *
 * The Stringify class provides the \ref format method. Its prototype is the same
 * as the printf function. As a result, it gives a string object holding the same
 * text that the printf function would have produced.
 *
 * \code
 *  cout << Stringify::format ("we got %d items which represents %f percent", 123, 57.3) << endl;
 * \endcode
 */
class Stringify
{
public:

    /** Generates a string object with a printf-like prototype.
     * \param[in] fmt : the format of the string
     * \param[in] ... : variable number of arguments
     * \return the string
     */
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
