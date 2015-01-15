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

#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
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
OptionsParser::OptionsParser (const string& name)
{
    setName (name);
    _proceed    = 0;
    _currentArg = 0;
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
void OptionsParser::add (Option* option)
{
    DEBUG (("OptionsParser::add  this=%p  option=%s\n", this, option->getLabel().c_str()));

    if (option != 0)
    {
        option->use ();

        /* We add the option in the list. */
        _options.push_back (option);
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
void OptionsParser::add (IOptionsParser* parser)
{
    /** UGGLY... */
    OptionsParser* p1 = dynamic_cast<OptionsParser*> (parser);
    if (p1)
    {
        for (list<Option*>::const_iterator it = p1->getOptions().begin(); it != p1->getOptions().end(); ++it)
        {
            this->add (*it);
        }
    }

    /** UGGLY... */
    OptionsParserComposite* p2 = dynamic_cast<OptionsParserComposite*> (parser);
    if (p2)
    {
        for (vector<IOptionsParser*>::const_iterator it = p2->getParsers().begin(); it != p2->getParsers().end(); ++it)
        {
            this->add (*it);
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
int OptionsParser::remove (const char * label)
{
    for (std::list<Option*>::iterator it = _options.begin(); it != _options.end(); it++)
    {
		if(  ! strcmp ((*it)->getLabel().c_str(), label) )
		{
			_options.erase(it);
			break;
		}
		
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

    _argIdxOk.resize (argc);

    /** We have to reset all options to their default values. */
    for (list<Option*>::iterator it = _options.begin();  it != _options.end();  it++)
    {
        (*it)->_param =  (*it)->_defaultParam;
    }

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
                    _warnings.push_back (make_pair(_currentArg-1,buffer));
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

                    /** We tag the arguments that have been correctly parsed. */
                    size_t startIdx = _currentArg-optionMatch->getNbArgs()-1;
                    for (size_t k=0; k<=optionMatch->getNbArgs(); k++)  {  _argIdxOk[startIdx+k] = true;  }

                    res=0;  // reduce warning
                }
            }
        }
        else
        {
            /** Unknown option => add it to the warning list. */
            _warnings.push_back (make_pair (_currentArg-1,Stringify::format("Unknown '%s'", txt)));
        }
    }

    /** We check mandatory options. */
    checkMandatoryOptions ();

    /* Now, we check that the accepted options are used with with include options. */
    checkIncludingOptions ();

    DEBUG (("OptionsParser::parse  errorsNb=%d  warningsNb=%d \n", _errors.size(), _warnings.size()));

    /** We may launch an exception if needed. */
    if (!_errors.empty())   {  throw OptionFailure (*this);  }

    /** We fill the properties. */
    buildProperties ();

	//if (!_warnings.empty())   {  throw OptionFailure (*this);  }

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
void OptionsParser::buildProperties ()
{
    setProperties (new Properties());

    for (std::list<Option*>::iterator it = _options.begin();  it != _options.end();  it++)
    {
        /** Shortcut. */
        Option* opt = *it;

        if (saw(opt->getLabel()) || opt->getParam().empty()==false)
        {
            _properties->add (0, (*it)->getLabel(), (*it)->getParam());
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
        snprintf (buffer, sizeof(buffer), "Too few arguments for the %s option...", option->getLabel().c_str());
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
void OptionsParser::displayErrors (std::ostream& os, size_t level)  const
{
    for (list<string>::const_iterator it = _errors.begin();  it != _errors.end();  it++)
    {
        os << /*indent(os,level)*/  "Error : " << it->c_str() << endl;
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
void OptionsParser::displayWarnings (std::ostream& os, size_t level, std::vector<bool>* idx)  const
{
    for (list<pair<size_t,string> >::const_iterator it = _warnings.begin();  it != _warnings.end();  it++)
    {
        bool shouldPrint = true;

        if (idx)
        {
            if ((*idx)[it->first]==false)  { shouldPrint = true;   (*idx) [it->first] = true;  }
            else                           { shouldPrint = false; }
        }

        if (shouldPrint)
        {
            os << /*indent(os,level)*/ "Warning : " << it->second.c_str() << endl;
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
void OptionsParser::displayHelp (std::ostream& os, size_t level)  const
{
    /** We first look for the longest option name. */
    size_t maxLen=0;
    for (list<Option*>::const_iterator it = _options.begin();  it != _options.end();  it++)
    {
        if (!(*it)->getLabel().empty())  {  maxLen = std::max (maxLen, (*it)->getLabel().size());  }
    }

    os << endl;
    indent(os,level) <<  "[" << _name << " options]" << endl;

    for (list<Option*>::const_iterator it = _options.begin();  it != _options.end();  it++)
    {
        if (!(*it)->getLabel().empty() && (*it)->isVisible())
        {
            if ((*it)->getNbArgs() > 0)
            {
                indent(os,level) << Stringify::format ("    %-*s (%d arg) :    %s",
                    (int)maxLen,
                    (*it)->getLabel().c_str(),
                    (int)(*it)->getNbArgs(),
                    (*it)->getHelp().c_str(),
                    (*it)->getDefaultParam().c_str()
                );

                if ((*it)->isMandatory()==false)  {  os << Stringify::format ("  [default '%s']", (*it)->getDefaultParam().c_str());  }

                os << endl;
            }
            else
            {
                indent(os,level) << Stringify::format ("    %-*s (%d arg) :    %s\n",
                    (int)maxLen,
                    (*it)->getLabel().c_str(),
                    (int)(*it)->getNbArgs(),
                    (*it)->getHelp().c_str()
                );
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

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParser::hide (const char* label)
{
	Option* opt = lookForOption ((char*)label);

	if (opt != 0)  { opt->hide (); }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
OptionsParserAbstract::OptionsParserAbstract () : _properties(0), _visible(true)
{
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
misc::IProperties* OptionsParserAbstract::parseString (const std::string& s)
{
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
    _properties = this->parse (argc, argv);

    DEBUG (("OptionsParser::parseString  argc=%d => idx=%ld\n", argc, idx));

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
void OptionsParserAbstract::displayVersion (std::ostream& os, size_t level)  const
{
    indent(os,level)  << Stringify::format ("* version %s (%s)\n* built on %s with compiler '%s'\n* supported kmer sizes %d %d %d %d",
        System::info().getVersion().c_str(),
        System::info().getBuildDate().c_str(),
        System::info().getBuildSystem().c_str(),
        System::info().getBuildCompiler().c_str(),
        KSIZE_1, KSIZE_2, KSIZE_3, KSIZE_4
    ) << endl;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
OptionsParserComposite::~OptionsParserComposite ()
{
    for (size_t i=0; i<_parsers.size(); i++)  {  _parsers[i]->forget(); }
    _parsers.clear();
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
misc::IProperties* OptionsParserComposite::parse (int argc, char* argv[])
{
    setProperties (new Properties ());

    bool hasExceptions = false;

    /** We loop over each known parser. */
    for (size_t i=0; i<_parsers.size(); i++)
    {
        try
        {
            /** We parse the options with the current parser. */
            _parsers[i]->parse (argc, argv);
        }
        catch (OptionFailure& e)
        {
            hasExceptions = true;
        }
        catch (Exception& e)
        {
            hasExceptions = true;
        }

        /** We add the found options to the global one. */
        _properties->add (0, _parsers[i]->getProperties());
    }

    _argIdxOk.resize(argc);

    /** We merge information of correctly parsed arguments. */
    for (size_t i=0; i<_parsers.size(); i++)
    {
        for (size_t j=0; j<argc; j++)  {  _argIdxOk[j] = _argIdxOk[j] | _parsers[i]->getParsedArgIndexes()[j];  }
    }

    /** We count the number of arguments correctly parsed. */
    size_t nbParsed=0;   for (size_t j=0; j<argc; j++)  {  nbParsed += _argIdxOk[j] ? 1 : 0;  }

    bool shouldSendException = (argc == 1) || (hasExceptions == true) || (hasExceptions == false &&  nbParsed != argc-1);

    /** The parsing is correct if we have arguments or we had no exception and if the number of parsed arguments is ok. */
    if (shouldSendException)
    {
        /** We set the parsing mask for all parsers. */
        for (size_t i=0; i<_parsers.size(); i++) { _parsers[i]->setParsedArgIndexes (_argIdxOk); }

        /** We throw an exception. */
        throw OptionFailure(*this);
    }

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
void OptionsParserComposite::displayErrors (std::ostream& os, size_t level)  const
{
    /** We loop over each known parser. */
    for (size_t i=0; i<_parsers.size(); i++)  {  _parsers[i]->displayErrors (os, level+1);  }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParserComposite::displayWarnings (std::ostream& os, size_t level, std::vector<bool>* idx) const
{
    std::vector<bool> warningsIdx = idx==0 ? getParsedArgIndexes() : *idx;

    /** We loop over each known parser. */
    for (size_t i=0; i<_parsers.size(); i++)  {  _parsers[i]->displayWarnings (os, level+1, &warningsIdx);  }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void OptionsParserComposite::displayHelp (std::ostream& os, size_t level) const
{
    if (isVisible() == true)
    {
        os << endl;
        indent(os,level) <<  "[" << _name << " options]" << endl;

        /** We loop over each known parser. */
        for (size_t i=0; i<_parsers.size(); i++) { if (_parsers[i]->isVisible())  {  _parsers[i]->displayHelp (os, level+1);  }  }
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
IOptionsParser* OptionsParserComposite::getParser (const std::string& name)
{
    IOptionsParser* result = 0;

    for (size_t i=0; result==0 && i<_parsers.size(); i++)
    {
        if (_parsers[i]->getName() == name)  { result = _parsers[i]; }
    }

    return result;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
