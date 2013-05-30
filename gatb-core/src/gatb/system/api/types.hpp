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

#include <sys/types.h>

#endif /* _GATB_CORE_SYSTEM_TYPES_HPP_ */
