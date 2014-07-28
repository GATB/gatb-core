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

/** \file Enums.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Enumerations
 */

#ifndef _GATB_CORE_TOOLS_MISC_ENUMS_HPP_
#define _GATB_CORE_TOOLS_MISC_ENUMS_HPP_

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

enum BloomKind  {  BLOOM_NONE, BLOOM_BASIC, BLOOM_CACHE, BLOOM_DEFAULT };

static void parse (const std::string& s, BloomKind& kind)
{
         if (s == "none")     { kind = BLOOM_NONE;  }
    else if (s == "basic")    { kind = BLOOM_BASIC;  }
    else if (s == "cache")    { kind = BLOOM_CACHE; }
    else if (s == "default")  { kind = BLOOM_CACHE; }
    else   { throw system::Exception ("bad Bloom kind '%s'", s.c_str()); }
}

static const char* toString (BloomKind kind)
{
    switch (kind)
    {
        case BLOOM_NONE:      return "none";
        case BLOOM_BASIC:     return "basic";
        case BLOOM_CACHE:     return "cache";
        case BLOOM_DEFAULT:   return "cache";
        default:        throw system::Exception ("bad Bloom kind %d", kind);
    }
}

/********************************************************************************/

enum DebloomKind  { DEBLOOM_NONE, DEBLOOM_ORIGINAL, DEBLOOM_CASCADING, DEBLOOM_DEFAULT };

static void parse (const std::string& s, DebloomKind& kind)
{
         if (s == "none")       { kind = DEBLOOM_NONE;      }
    else if (s == "original")   { kind = DEBLOOM_ORIGINAL;  }
    else if (s == "cascading")  { kind = DEBLOOM_CASCADING; }
    else if (s == "default")    { kind = DEBLOOM_CASCADING; }
    else   { throw system::Exception ("bad debloom kind '%s'", s.c_str()); }
}

static std::string toString (DebloomKind kind)
{
    switch (kind)
    {
        case DEBLOOM_NONE:      return "none";
        case DEBLOOM_ORIGINAL:  return "original";
        case DEBLOOM_CASCADING: return "cascading";
        case DEBLOOM_DEFAULT:   return "cascading";
        default:        throw system::Exception ("bad debloom kind %d", kind);
    }
}

/********************************************************************************/

enum BranchingKind  { BRANCHING_NONE, BRANCHING_STORED };

static void parse (const std::string& s, BranchingKind& kind)
{
         if (s == "none")     { kind = BRANCHING_NONE;  }
    else if (s == "stored")   { kind = BRANCHING_STORED;  }
    else   { throw system::Exception ("bad branching kind '%s'", s.c_str()); }
}

static std::string toString (BranchingKind kind)
{
    switch (kind)
    {
        case BRANCHING_NONE:     return "none";
        case BRANCHING_STORED:   return "stored";
        default:        throw system::Exception ("bad branching kind %d", kind);
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_ENUMS_HPP_ */
