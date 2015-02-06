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

/** \file LibraryInfo.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_LIBRARY_INFO_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_LIBRARY_INFO_HPP_

/********************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Framework class for implementing tools (ie. binary tools).
 */
class LibraryInfo
{
public:

    /** Get information about the GATB-CORE library
     * \return information as a IProperties instance
     */
    static IProperties& getInfo()
    {
        static system::SmartObject singleton;

        if (singleton.hasRef() == false)
        {
            IProperties* props = new Properties();

            props->add (0, "gatb-core-library", "");
            props->add (1, "version",        "%s", system::impl::System::info().getVersion().c_str());
            props->add (1, "build_date",     "%s", system::impl::System::info().getBuildDate().c_str());
            props->add (1, "build_system",   "%s", system::impl::System::info().getBuildSystem().c_str());
            props->add (1, "build_compiler", "%s", system::impl::System::info().getBuildCompiler().c_str());
            //props->add (1, "build_options",  "%s", system::impl::System::info().getBuildOptions().c_str());
            props->add (1, "build_kmer_size", "%d %d %d %d", KSIZE_1, KSIZE_2, KSIZE_3, KSIZE_4);
            props->add (1, "custom_memalloc", "%d", CUSTOM_MEM_ALLOC);

            singleton.setRef (props);
        }
        return * (dynamic_cast<IProperties*>(singleton.getRef()));
    }

    /** Display information about the GATB-CORE library
     * \param[in] os : output stream used for dumping library information
     */
    static void displayVersion (std::ostream& os)
    {
        os << Stringify::format ("* version %s (%s)\n* built on %s with compiler '%s'\n* supported kmer sizes %d %d %d %d",
            system::impl::System::info().getVersion().c_str(),
            system::impl::System::info().getBuildDate().c_str(),
            system::impl::System::info().getBuildSystem().c_str(),
            system::impl::System::info().getBuildCompiler().c_str(),
            KSIZE_1, KSIZE_2, KSIZE_3, KSIZE_4
        ) << std::endl;
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_LIBRARY_INFO_HPP_ */
