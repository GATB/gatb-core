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

/** \file StringLine.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool framework
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_STRINGLINE_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_STRINGLINE_HPP_

/********************************************************************************/

#include <string>
#include <iostream>
#include <cstdarg>
#include <cstdio>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief String line formatter
 */
class StringLine
{
public:

    static std::string format (const std::string& str)
    {
        using namespace std;
        string spaces (getDefaultWidth(), ' ');
        string result;

        if (str.size() > spaces.size())  { result = str.substr (0,spaces.size()-getSplitWidth()) + string(getSplitWidth(),'.'); }
        else                             { result = (str + spaces).substr (0,getDefaultWidth());  }

        return result;
    }

    static std::string format (const char* fmt, ...)
    {
        using namespace std;

        /** We get the buffer from the ellipsis argument. */
        char* buffer = 0;
        va_list args;
        va_start (args, fmt);
        vasprintf (&buffer, fmt, args);
        va_end (args);

        string result (buffer != NULL ? buffer : "");

        if (buffer != NULL)  { free(buffer); }

        /** We return the result. */
        return format (result);
    }

    static size_t& getDefaultWidth () { static size_t value=40; return value; }
    static void setDefaultWidth (size_t value) { getDefaultWidth() = value; }

    static size_t& getSplitWidth () { static size_t value=3; return value; }
    static void setSplitWidth (size_t value) { getDefaultWidth() = value; }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_STRINGLINE_HPP_ */
