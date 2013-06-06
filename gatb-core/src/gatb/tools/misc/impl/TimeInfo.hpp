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

/** \file TimeTools.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for time measurement feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_TIMETOOLS_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_TIMETOOLS_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/system/api/ITime.hpp>

#include <map>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Tool for time statistics.
 *
 * This class provides methods for getting time information between two execution points.
 *
 * One can use a label for getting a specific time duration; it is possible later to get
 * this duration by giving the label.
 *
 * Example of use:
 * \code
 void foo ()
 {
     TimeInfo t;

     t.addEntry ("part1");
     // do something here
     t.stopEntry ("part1");

     t.addEntry ("part2");
     // do something here
     t.stopEntry ("part2");

     // now, we dump the duration of part1 and part2:
     cout << "part1: " << t.getEntryByKey("part1") << "  "
          << "part2: " << t.getEntryByKey("part2") << endl;
 }
 * \endcode
  */
class TimeInfo : public dp::SmartPointer
{
public:

    /** Default constructor. */
    TimeInfo ();

    /** Constructor taking a time factory.
     * \param[in] aTime : the time factory to be used.
     */
    TimeInfo (system::ITime& aTime);

    /** Get the start time for a given label.
     * \param[in] name : the label
     */
    virtual void start (const char* name);

    /** Get the stop time for a given label.
     * \param[in] name : the label
     */
    virtual void stop (const char* name);

    /** Provides (as a map) all got durations for each known label/
     * \return a map holding all retrieved timing information.
     */
    const std::map <std::string, u_int32_t>& getEntries ();

    /** Retrieve the duration for a given label.
     * \param[in] key : the label we want the duration for.
     * \return the duration.
     */
    u_int32_t getEntryByKey (const std::string& key);

    /** Creates and return as a IProperties instance the whole timing information.
     * \param[in] root : root name of the properties to be returned.
     * \return the created IProperties instance.
     */
    virtual tools::misc::IProperties* getProperties (const std::string& root);

private:

    system::ITime&  _time;
    std::map <std::string, u_int32_t>  _entriesT0;
    std::map <std::string, u_int32_t>  _entries;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_TIMETOOLS_HPP_ */
