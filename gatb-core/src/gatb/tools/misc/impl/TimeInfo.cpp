/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/system/impl/System.hpp>

#define DEBUG(a)  //printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
TimeInfo::TimeInfo () : _time(system::impl::System::time())
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
TimeInfo::TimeInfo (system::ITime& aTime) : _time(aTime)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void TimeInfo::start (const char* name)
{
    _entriesT0 [name] = _time.getTimeStamp();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void TimeInfo::stop (const char* name)
{
    _entries [name] += _time.getTimeStamp() - _entriesT0 [name];
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
const std::map <std::string, u_int32_t>& TimeInfo::getEntries ()
{
    return _entries;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
u_int32_t TimeInfo::getEntryByKey (const std::string& key)
{
    u_int32_t result = 0;
    std::map <std::string, u_int32_t>::iterator  it = _entries.find (key);
    if (it != _entries.end())  { result = it->second; }
    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
tools::misc::IProperties* TimeInfo::getProperties (const std::string& root)
{
    tools::misc::IProperties* props = new tools::misc::impl::Properties();

    props->add (0, root,   "");

    std::map <std::string, u_int32_t>::const_iterator  it;
    for (it = getEntries().begin(); it != getEntries().end();  it++)
    {
        props->add (1, it->first.c_str(), "%ld", it->second);
    }

    return props;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
