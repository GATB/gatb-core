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

#include <cstdarg>
#include <cstdio>

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
struct HierarchyParserVisitor : public IOptionsParserVisitor
{
    void visitOptionsParser (impl::OptionsParser& object, size_t depth)
    {
        for (std::list<IOptionsParser*>::const_iterator it = object.getParsers().begin(); it != object.getParsers().end(); ++it)
        {
            (*it)->accept (*this, depth+1);
        }
    }
};

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
struct PostParserVisitor : public HierarchyParserVisitor
{
    PostParserVisitor (const set<string>& foundParsers, IOptionsParser::Result& result)
        : foundParsers(foundParsers), result(result) {}

    void visitOption (Option& object, size_t depth)
    {
        /** The current Option may have been not seen during the parsing. */
        if (foundParsers.find(object.getName()) == foundParsers.end())
        {
            /** This option is mandatory, so this is an error. */
            if (object.isMandatory()==true)
            {
                result.errors.push_back (Stringify::format ("Option '%s' is mandatory", object.getName().c_str()));
            }
            /** This option is not mandatory, so we use its default value. */
            else if (object.getDefaultValue().empty()==false)
            {
                result.properties.add (0, object.getName(), object.getDefaultValue());
            }
        }
    }

    const set<string>&      foundParsers;
    IOptionsParser::Result& result;
};

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
struct ParserVisitor : public IOptionsParserVisitor
{
    int argc; char** argv;
    IOptionsParser::Result result;

    ParserVisitor (int argc, char** argv) : argc(argc), argv(argv) {}

    IOptionsParser::Result& getResult () { return result; }

    void visitOptionsParser (OptionsParser& object, size_t depth)
    {
        set<string> foundParsers;

        /** We loop the arguments. */
        for (argc--,argv++; argc>0; argc--, argv++)
        {
            char* txt = argv[0];

            IOptionsParser* match = object.getParser (txt);

            if (match != 0)
            {
                /** We parse recursively with the found parser. */
                match->accept (*this, depth+1);

                /** We keep a track of found parsers. */
                foundParsers.insert (match->getName());
            }
            else
            {
                result.errors.push_back (Stringify::format("Unknown '%s'", txt));
            }
        }

        PostParserVisitor visitor (foundParsers, result);
        object.accept (visitor);
    }

    void visitOption (Option& object, size_t depth)
    {
        /** We go to the first argument of the current option. */
        argc--; argv++;

        /** We retrieve the arguments of the option. */
        std::list<std::string> optionArgs;
        for (int i=0; argc>0; argc--, argv++)
        {
            optionArgs.push_back (argv[0]);

            if (++i>=object.getNbArgs())  { break; }
        }

        if (optionArgs.size() == object.getNbArgs())
        {
            object.proceed (optionArgs, result.properties);
        }
        else
        {
            char buffer [128];
            snprintf (buffer, sizeof(buffer), "Too few arguments for the %s option...", object.getName().c_str());
            result.errors.push_back (buffer);
        }
    }
};


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
OptionsParser::OptionsParser (const std::string& name)
    : _name(name), _visible(true), _properties(0)
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
OptionsParser::~OptionsParser ()
{
    for (list<IOptionsParser*>::const_iterator it = _parsers.begin(); it != _parsers.end(); ++it)
    {
        (*it)->forget();
    }

    setProperties (0);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
misc::IProperties* OptionsParser::parse (int argc, char** argv)
{
    /** We parse the arguments through a visitor. */
    ParserVisitor visitor (argc, argv);
    this->accept (visitor);

    /** We launch an exception in case we got errors during the parsing. */
    if (!visitor.getResult().errors.empty())  {  throw OptionFailure (this, visitor.getResult());  }

    /** We set the properties. */
    setProperties (new Properties (visitor.getResult().properties));

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
misc::IProperties* OptionsParser::parseString (const std::string& s)
{
    int    argc   = 0;
    size_t idx    = 0;

    /** We tokenize the string and count the number of tokens. */
    TokenizerIterator it1 (s.c_str(), " ");
    for (it1.first(); !it1.isDone(); it1.next())  { argc++; }

    argc++;

    /** We allocate a table of char* */
    char** argv = (char**) calloc (argc, sizeof(char*));

    argv[idx++] = strdup (getName().c_str());

    /** We fill the table with the tokens. */
    TokenizerIterator it2 (s.c_str(), " ");
    for (it2.first(); !it2.isDone(); it2.next())
    {
        argv[idx++] = strdup (it2.item());
    }

    for (int i=0; i<argc; i++)  {  DEBUG (("   item='%s' \n", argv[i]));  }

    /** We parse the arguments. */
    misc::IProperties* result = this->parse (argc, argv);

    DEBUG (("OptionsParser::parseString  argc=%d => idx=%ld\n", argc, idx));

    /** Some clean up. */
    for (int i=0; i<argc; i++)  {  free (argv[i]);  }
    free (argv);

    /** We return the result. */
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
bool OptionsParser::saw (const std::string& name) const
{
    return (_properties != 0  && _properties->get(name) != 0);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
IOptionsParser* OptionsParser::getParser (const std::string& name)
{
    DEBUG (("OptionsParser::getParser  '%s'  look for '%s'  parsers.size=%ld   equal=%d\n",
        getName().c_str(), name.c_str(), _parsers.size(), name == getName()
    ));

    if (name == getName()) { return  this; }

    IOptionsParser* result = 0;
    for (list<IOptionsParser*>::iterator it = _parsers.begin(); result==0 && it != _parsers.end(); ++it)
    {
        result = (*it)->getParser(name);
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
class OptionsHelpVisitor : public HierarchyParserVisitor
{
public:

    OptionsHelpVisitor (std::ostream& os) : os(os),nameMaxLen(0) {}

    void visitOptionsParser (OptionsParser& object, size_t depth)
    {
        if (object.isVisible() == true)
        {
            /** We first look for the longest option name. */
            nameMaxLen=0;
            for (list<IOptionsParser*>::const_iterator it = object.getParsers().begin(); it != object.getParsers().end(); ++it)
            {
                if (!(*it)->getName().empty())  {  nameMaxLen = std::max (nameMaxLen, (*it)->getName().size());  }
            }

            os << endl;
            indent(os,depth) <<  "[" << object.getName() << " options]" << endl;

            /** We loop over each known parser. */
            for (list<IOptionsParser*>::const_iterator it = object.getParsers().begin(); it != object.getParsers().end(); ++it)
            {
                if ((*it)->isVisible())  {  (*it)->accept (*this, depth+1); }
            }
        }
    }

    void visitOption (Option& object, size_t depth)
    {
        if (!object.getName().empty() && object.isVisible())
        {
            if (object.getNbArgs() > 0)
            {
                indent(os,depth) << Stringify::format ("    %-*s (%d arg) :    %s",
                    (int)nameMaxLen,
                    object.getName().c_str(),
                    (int)object.getNbArgs(),
                    object.getHelp().c_str(),
                    object.getDefaultValue().c_str()
                );

                if (object.isMandatory()==false)  {  os << Stringify::format ("  [default '%s']", object.getDefaultValue().c_str());  }

                os << endl;
            }
            else
            {
                indent(os,depth) << Stringify::format ("    %-*s (%d arg) :    %s\n",
                    (int)nameMaxLen,
                    object.getName().c_str(),
                    (int)object.getNbArgs(),
                    object.getHelp().c_str()
                );
            }
        }
    }

private:

    std::ostream& os;
    size_t        nameMaxLen;

    std::ostream& indent (std::ostream& os, size_t level)  const { for (size_t i=0; i<level; i++)  { os << "   "; }  return os; }
};

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int OptionFailure::displayErrors (std::ostream& os) const
{
    for (std::list<std::string>::const_iterator it = _result.errors.begin(); it != _result.errors.end(); ++it)
    {
        os << "ERROR: " << *it << std::endl;
    }

    OptionsHelpVisitor visitor (os);
    _parser->accept (visitor, 0);

    return EXIT_FAILURE;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
