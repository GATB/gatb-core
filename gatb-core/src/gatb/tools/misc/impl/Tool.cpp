/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#define DEBUG(a)  //printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

const char* Tool::STR_NB_CORES   = "-nb-cores";
const char* Tool::STR_QUIET      = "-quiet";
const char* Tool::STR_STATS_XML  = "-stats";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Tool::Tool (const std::string& name) : _name(name), _input(0), _output(0), _parser(0), _dispatcher(0)
{
    setOutput (new Properties());
    _output->add (0, _name.c_str());

    setParser (new OptionsParser ());
    _parser->add (new OptionOneParam (Tool::STR_NB_CORES,   "number of cores",                      false));
    _parser->add (new OptionNoParam  (Tool::STR_QUIET,      "quiet execution",                      false));
    _parser->add (new OptionOneParam (Tool::STR_STATS_XML,  "dump exec info into a XML file",       false));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Tool::~Tool ()
{
    setInput      (0);
    setOutput     (0);
    setParser     (0);
    setDispatcher (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperties* Tool::run (int argc, char* argv[])
{
    return run (getOptionsParser()->parse (argc, argv));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperties* Tool::run (IProperties* input)
{
    /** We keep the input parameters. */
    setInput (input);

    /** We add a potential config file to the input properties. */
    _input->add (1, new Properties (System::info().getHomeDirectory() + "/." + getName() ));

    /** We add the input properties to the output result. */
    _output->add (1, input);

    /** We define one dispatcher. */
    setDispatcher (new ParallelCommandDispatcher (input->getInt(STR_NB_CORES)));

    /** We execute the actual job. */
    {
        //TIME_INFO (_timeInfo, _name);
        execute ();
    }

    /** We add the time properties to the output result. */
    _output->add (1, _timeInfo.getProperties ("time"));

    /** We may have some post processing. */
    postExecute ();

    /** We return the output properties. */
    return _output;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Tool::postExecute ()
{
    /** We may have to dump execution information into a stats file. */
    if (_input->get(Tool::STR_STATS_XML) != 0)
    {
        XmlDumpPropertiesVisitor visit (_output->getStr (Tool::STR_STATS_XML));
        _output->accept (&visit);
    }

    /** We may have to dump execution information to stdout. */
    if (_input->get(Tool::STR_QUIET) == 0)
    {
        RawDumpPropertiesVisitor visit;
        _output->accept (&visit);
    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
