/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#define DEBUG(a)  printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

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
Algorithm::Algorithm (const std::string& name, gatb::core::tools::misc::IProperties* input)
    : _name(name), _input(0), _output(0), _info(0), _dispatcher(0)
{
    setInput  (input ? input : new Properties());
    setOutput (new Properties());
    setInfo   (new Properties());
    setDispatcher (new ParallelCommandDispatcher (_input->get(STR_NB_CORES)  ? _input->getInt(STR_NB_CORES) : 0) );

    _info->add (0, _name);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Algorithm::~Algorithm ()
{
    setInput      (0);
    setOutput     (0);
    setInfo       (0);
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
dp::IteratorListener* Algorithm::createIteratorListener (size_t nbIterations, const char* message)
{
    if (getInput()->get(STR_PROGRESS_BAR)==0)  { return new IteratorListener(); }

    switch (getInput()->getInt(STR_PROGRESS_BAR))
    {
        case 0: default:    return new IteratorListener ();
        case 1:             return new Progress      (nbIterations, message);
        case 2:             return new ProgressTimer (nbIterations, message);
    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
