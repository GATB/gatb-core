/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>
#include <gatb/system/impl/System.hpp>

#include <stdarg.h>
#include <stdio.h>

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
OptionsParser::OptionsParser ()  : _properties(0)
{
    _proceed=0;
    setProperties (new Properties());
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
OptionsParser::~OptionsParser ()
{
    setProperties (0);

    for (std::list<Option*>::iterator it = _options.begin(); it != _options.end(); it++)
    {
        (*it)->forget();
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int OptionsParser::add (Option* option)
{
    DEBUG (("OptionsParser::add  this=%p  option=%s\n", this, option->getLabel().c_str()));

    if (option != 0)
    {
        option->use ();

        /* We add the option in the list. */
        _options.push_back (option);
    }

    return _options.size();
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
misc::IProperties* OptionsParser::parse (int argc, char* argv[])
{
    DEBUG (("OptionsParser::parse  this=%p  argc=%d\n", this, argc));

    /* Some initializations. */
    _argc=argc; _argv=argv; _currentArg=1;
    _proceed=1;
    _errors.clear();
    _warnings.clear();
    _seenOptions.clear();
    char* exclude;

    char* txt=0;
    while ( (txt = nextArg ()) != 0)
    {
        DEBUG (("OptionsParser::parse:  _currentArg=%d  txt=%p \n", _currentArg, txt));

        /* We look if is matches one of the recognized options. */
        Option* optionMatch = lookForOption (txt);
        DEBUG (("OptionsParser::parse:  txt='%s'  optionMatch=%p \n", txt, optionMatch));

        if (optionMatch != 0)
        {
            /* This is a recognized option. So we try to retrieve the args that should follow. */
            std::list<std::string> optionArgs;
            getOptionArgs (optionMatch, optionArgs);

            DEBUG (("OptionsParser::parse:  optionMatch.size=%ld  optionArgs.size=%ld\n",
                optionMatch->getNbArgs(),
                optionArgs.size()
            ));

            if (optionArgs.size() == optionMatch->getNbArgs())
            {
                /* We first look if this option has not already been met. */
                if (!optionMatch->canBeMultiple () && saw (optionMatch->getLabel()))
                {
                    char buffer[128];
                    snprintf (buffer, sizeof(buffer), "Option %s already seen, ignored...", optionMatch->getLabel().c_str());
                    _warnings.push_back (buffer);
                }
                else if ( (exclude = checkExcludingOptions (optionMatch)) != 0)
                {
                    char buffer[128];
                    snprintf (buffer, sizeof(buffer), "Option %s can't be used with option %s, ignored...",
                        optionMatch->getLabel().c_str(),
                        exclude
                    );
                    _errors.push_back (buffer);
                }
                else
                {
                    int res = optionMatch->proceed (optionArgs);
                    _seenOptions.push_back (optionMatch);

                    DEBUG (("OptionsParser::parse:  proceed the option => res=%ld  seenOptions=%ld  _currentArg=%d\n",
                        res, _seenOptions.size(), _currentArg
                    ));

                    res=0;  // reduce warning
                }
            }
        }
        else
        {
#if 0
            /* This is an NOT a recognized option. We try to find an No_Option. */
            giveToNoOption (txt);
#else
            /** Unknown option => add it to the warning list. */
            _warnings.push_back (txt);
#endif
        }
    }

    /** We check mandatory options. */
    checkMandatoryOptions ();

    /* Now, we check that the accepted options are used with with include options. */
    checkIncludingOptions ();

    DEBUG (("OptionsParser::parse  errorsNb=%d  warningsNb=%d \n", _errors.size(), _warnings.size()));

    /** We may launch an exception if needed. */
    if (!_errors.empty())   {  throw OptionFailure (this);  }

    /** We fill the properties. */
    buildProperties ();

    /* And we return the errors number. */
    return _properties;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
misc::IProperties* OptionsParser::parse (const std::string& s)
{
    IProperties* result = 0;
    int    argc   = 0;
    size_t idx    = 0;

    /** We tokenize the string and count the number of tokens. */
    TokenizerIterator it1 (s.c_str(), " ");
    for (it1.first(); !it1.isDone(); it1.next())  { argc++; }

    argc++;

    /** We allocate a table of char* */
    char** argv = (char**) calloc (argc, sizeof(char*));

    argv[idx++] = strdup ("foo");

    /** We fill the table with the tokens. */
    TokenizerIterator it2 (s.c_str(), " ");
    for (it2.first(); !it2.isDone(); it2.next())
    {
        argv[idx++] = strdup (it2.item());
    }

    for (int i=0; i<argc; i++)  {  DEBUG (("   item='%s' \n", argv[i]));  }

    /** We parse the arguments. */
    result = parse (argc, argv);

    DEBUG (("OptionsParser::parse  argc=%d => idx=%ld  result=%d  seenOptions=%ld\n", argc, idx, result, _seenOptions.size() ));

    /** Some clean up. */
    for (int i=0; i<argc; i++)  {  free (argv[i]);  }
    free (argv);

    /** We return the result. */
    return _properties;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::buildProperties ()
{
    _properties->add (0, "input");

    for (std::list<Option*>::iterator it = _options.begin();  it != _options.end();  it++)
    {
        /** Shortcut. */
        Option* opt = *it;

        if (saw(opt->getLabel()) || opt->getParam().empty()==false)
        {
            _properties->add (1, (*it)->getLabel(), (*it)->getParam());
        }
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
char* OptionsParser::nextArg ()
{
    DEBUG (("\nCheckOption::nextArg:  _currentArg=%d  _argc=%d \n", _currentArg, _argc));

    return (_currentArg >= _argc ? (char*)0 : _argv[_currentArg++]);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
Option* OptionsParser::lookForOption (char* txt)
{
    Option* result = 0;

    DEBUG (("CheckOption::lookForOption:  txt='%s'  _options.size=%ld  \n", txt, _options.size()));

    for (list<Option*>::iterator it = _options.begin(); !result &&  it != _options.end(); it++)
    {
        if ((*it)->getLabel().compare (txt) == 0)  {  result = *it; }
    }

    DEBUG (("CheckOption::lookForOption:  _options.size=%ld  => found=%d \n", _options.size(), result!=0));

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
void OptionsParser::getOptionArgs (const Option* option, std::list<std::string>& result)
{
    char* txt;
    int i=1;
    int n = option->getNbArgs();
    while ( (i<=n) && (txt=nextArg()) )
    {
        result.push_back (txt);
        i++;
    }

    if (i<=n)
    {
        char buffer [128];
        snprintf (buffer, sizeof(buffer), "To few arguments for the %s option...", option->getLabel().c_str());
        _errors.push_back (buffer);
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::displayErrors (FILE* fp)
{
    for (list<string>::iterator it = _errors.begin();  it != _errors.end();  it++)
    {
        fprintf (fp, "Error : %s\n", it->c_str());
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::displayWarnings (FILE* fp)
{
    for (list<string>::iterator it = _warnings.begin();  it != _warnings.end();  it++)
    {
        fprintf (fp, "Warning : %s\n", it->c_str());
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::displayHelp (FILE* fp)
{
    for (list<Option*>::iterator it = _options.begin();  it != _options.end();  it++)
    {
        if (!(*it)->getLabel().empty())
        {
            fprintf (fp, "\t%s :\t %s\n",  (*it)->getLabel().c_str(),  (*it)->getHelp().c_str());
        }
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::giveToNoOption (char* txt)
{
    Option* option = 0;

    for (list<Option*>::iterator it = _options.begin();  option==0  &&  it != _options.end();  it++)
    {
        if ((*it)->getLabel().empty())  {  option = (*it);  }
    }

    if (option != 0)
    {
        list<string> args;
        args.push_back (txt);
        option->proceed (args);
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
bool OptionsParser::saw (const std::string& txt)
{
    bool found = false;

    for (std::list<Option*>::iterator it = _seenOptions.begin();  !found && it != _seenOptions.end();  it++)
    {
        found = (*it)->getLabel().compare(txt) == 0;
    }

    DEBUG (("CheckOption::saw: txt='%s'  _seenOptions.size=%ld  => found=%d \n",
        txt.c_str(), _seenOptions.size(), found
    ));

    return found;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
char* OptionsParser::checkExcludingOptions (const Option* option)
{
    if (option->getExclude().empty())
    {
        return (char*)0;
    }

    short found=0;

    for (list<Option*>::iterator it = _seenOptions.begin(); it != _seenOptions.end(); it++)
    {
        if (strstr ((*it)->getExclude().c_str(), option->getLabel().c_str()))
        {
            found=1;
            break;
        }
    }
    return (found ? (char*) option->getLabel().c_str() : (char*)0);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::checkIncludingOptions ()
{
    for (list<Option*>::iterator it = _seenOptions.begin(); it != _seenOptions.end(); it++)
    {
        string include = (*it)->getInclude ();
        if (! include.empty ())
        {
            short inner_found=0;

            for (list<Option*>::iterator itInner = _seenOptions.begin(); itInner != _seenOptions.end(); itInner++)
            {
                if (strstr (include.c_str(), (*itInner)->getLabel().c_str()))
                {
                    inner_found=1;
                    break;
                }
            }
            if (!inner_found)
            {
                char buffer[128];
                snprintf (buffer, sizeof(buffer),
                    "Option %s must be used with one of the following options : %s",
                    (*it)->getLabel().c_str(),
                    include.c_str()
                );

                _errors.push_back (buffer);
            }
        }
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::checkMandatoryOptions ()
{
    for (list<Option*>::iterator it = _options.begin(); it != _options.end(); it++)
    {
        bool found = false;

        for (list<Option*>::iterator itInner = _seenOptions.begin(); !found &&  itInner != _seenOptions.end(); itInner++)
        {
            found = (*it)->getLabel() == (*itInner)->getLabel();
        }

        DEBUG (("OptionsParser::checkMandatoryOptions : label=%s  mandatory=%d  found=%d\n",
            (*it)->getLabel().c_str(),
            (*it)->isMandatory(),
            found
        ));

        if ((*it)->isMandatory() && !found)
        {
            char buffer[128];
            snprintf (buffer, sizeof(buffer), "Option '%s' is mandatory", (*it)->getLabel().c_str());
            _errors.push_back (buffer);
        }
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
const Option* OptionsParser::getSeenOption (const std::string& label)
{
    const Option* result = 0;

    for (list<Option*>::iterator it = _seenOptions.begin();  !result && it != _seenOptions.end(); it++)
    {
        if ((*it)->getLabel().compare (label) == 0)   { result = (*it); }
    }

    return result;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
