/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Tool.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool framework
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Framework class for implementing tools (ie. binary tools).
 */
class Tool
{
public:

    /** Constructor.
     * \param[in] name: name of the tool. */
    Tool (const std::string& name);

    /** Destructor. */
    virtual ~Tool ();

    /** Get tool name
     * \return the tool name. */
    std::string getName () const  { return _name; }

    /** Run the tool with input parameters provided as a IProperties instance
     * \param[in] input : input parameters
     * \return
     */
    IProperties* run (IProperties* input);

    /** */
    IProperties* run (int argc, char* argv[]);

    /** */
    virtual OptionsParser* getOptionsParser ()  {  return _parser;  }

    /** */
    static const char* STR_NB_CORES;
    static const char* STR_QUIET;
    static const char* STR_STATS_XML;

protected:

    /** */
    virtual void execute () = 0;

    /** */
    virtual void postExecute ();

    std::string _name;

    IProperties* _input;
    void setInput (IProperties* input)  { SP_SETATTR(input); }

    IProperties* _output;
    void setOutput (IProperties* output)  { SP_SETATTR(output); }

    OptionsParser* _parser;
    void setParser (OptionsParser* parser)  { SP_SETATTR(parser); }

    dp::ICommandDispatcher* _dispatcher;
    void setDispatcher (dp::ICommandDispatcher* dispatcher)  { SP_SETATTR(dispatcher); }

    TimeInfo _timeInfo;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_ */
