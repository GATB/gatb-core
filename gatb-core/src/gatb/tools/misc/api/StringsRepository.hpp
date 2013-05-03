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

/** \file StringsRepository.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Pool of strings providing information to end users
 */

#ifndef _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_
#define _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Tools package */
namespace tools     {
/** \brief Misc interfaces */
namespace misc      {
/********************************************************************************/

/** \brief Pool of strings
 *
 * This class provides constant strings used throughout the code. It may be interesting
 * in order to have a central point for constant strings management:
 *      - ease translation in different languages
 *      - entry point for strings obfuscation if needed
 *
 * It could also be possible to read the strings from a configuration file.
 *
 * The class defines one (static) method per constant string to be used. Note that we
 * also (see below) define a macro definition that eases the use of such a facility.
 */

class StringRepository
{
public:

    /** \brief Singleton method.
     *
     * This method could return different types for the string repository, for translation
     * for instance.
     *
     * \return the singleton instance.
     */
    static StringRepository& singleton()  { static StringRepository instance; return instance; }

    const char* BANK_bad_file_number    () { return "bank files number is %d but should be in [1..%d]"; }
    const char* BANK_bad_file_path      () { return "unable to find file '%s'"; }
    const char* BANK_unable_open_file   () { return "error opening file: %s"; }
    const char* BANK_unable_write_file  () { return "unable to write into file"; }

};

/********************************************************************************/

#define STR_BANK_bad_file_number    gatb::core::tools::misc::StringRepository::singleton().BANK_bad_file_number ()
#define STR_BANK_bad_file_path      gatb::core::tools::misc::StringRepository::singleton().BANK_bad_file_path ()
#define STR_BANK_unable_open_file   gatb::core::tools::misc::StringRepository::singleton().BANK_unable_open_file ()
#define STR_BANK_unable_write_file  gatb::core::tools::misc::StringRepository::singleton().BANK_unable_write_file ()

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_ */
