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

/** \file OptionsParser.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for parsing command line arguments
 */

#ifndef _GATB_CORE_TOOLS_MISC_OPTION_PARSER_HPP_
#define _GATB_CORE_TOOLS_MISC_OPTION_PARSER_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>

#include <list>
#include <string>
#include <stdio.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Interface for command line option
 *
 * Define an interface for what is and what contains an option.
 *
 * It can't be instantiated since the constructor is protected and so, has to be derived.
 *
 * \see OptionsParser
 */
class Option : public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] name : name of the option
     * \param[in] nbArgs : number of arguments for this option
     * \param[in] mandatory : tells whether this option is mandatory or not
     * \param[in] defaultValue : default value for the option
     * \param[in] help : textual help for this option
     * \param[in] multiple : tells whether this option may be used more than once
     * \param[in] include : list of names of options that must be used with the current one
     * \param[in] exclude : list of names of options that must not be used with the current one
     */
    Option (
        const std::string& name,
        int nbArgs,
        bool mandatory,
        const std::string& defaultValue,
        const std::string& help,
        int multiple,
        const std::string& include,
        const std::string& exclude
    )
        : _name(name), _nbArgs(nbArgs), _mandatory(mandatory), _help(help), _multiple(multiple), _include(include), _exclude(exclude), _param(defaultValue)
    {
    }

    /** Desctructor. */
    virtual ~Option() {}

    /** Label of the option (example "-log").
     * \return the label
     */
    const std::string& getLabel () const   { return _name; }

    /** Parameter string
     * \return the parameter.
     */
    const std::string& getParam ()  const { return _param; }

protected:

    /** Gives the number of arguments that must follow the option.
     * \return the arguments number.
     */
    size_t getNbArgs () const   { return _nbArgs; }

    /** Tells whether the option is mandatory or not.
     * \return the mandatory status.
     */
    bool isMandatory () const { return _mandatory; }

    /** Help text about this option.
     * \return help string
     */
    const std::string& getHelp () const { return _help; }

    /** Tell if the option can be used more than once.
     * \return true if can be multiple, false otherwise
     */
    short canBeMultiple () const    { return _multiple; }

    /** List of options that should be used with the current one.
     * The format of this list is for example "-C,-x,-X")
     * \return list of options
     */
    const std::string& getInclude () const  { return _include; }

    /** List of options that must not be used with the current one.
     * The format of this list is for example "-F,-x,-X")
     * \return list of options
     */
    const std::string& getExclude ()  const { return _exclude; }

    /** When an option is recognized in the argumenst list, we look the number of waited args and put
     * them in a list of string objects. This is this list that is given as argument of the proceed() method
     * that mainly will affect the given args to the variable given to the instanciation of the
     * (derived class) Option.
     */
    virtual int proceed (const std::list<std::string>& args) = 0;

    std::string     _name;
    size_t          _nbArgs;
    bool            _mandatory;
    std::string     _help;
    short           _multiple;
    std::string     _include;
    std::string     _exclude;
    std::string     _param;

    /* Since the CheckOption class is responsable to the full job, we let it access to the internal informations
     of one Option. */
    friend class OptionsParser;
};

/********************************************************************************/

/** \brief Option that has no argument.
 *
 * This is a special option (with no name) that memorize the arguments that are not
 * involved with a known option.
 */
class OptionNoParam : public Option
{
public:

    /** Constructor.
     * \param[in] name : name of the option
     * \param[in] help : textual help for this option
     * \param[in] mandatory : tells whether this option is mandatory or not
     * \param[in] multiple : tells whether this option may be used more than once
     * \param[in] include : list of names of options that must be used with the current one
     * \param[in] exclude : list of names of options that must not be used with the current one
     */
    OptionNoParam (
        const std::string& name,
        const std::string& help,
        bool mandatory = false,
        int multiple = 0,
        const char* include = "",
        const char* exclude = ""
    )
        : Option (name, 0, mandatory, "", help, multiple, include, exclude)
    {
    }

    /** \copydoc Option::proceed */
    int proceed (const std::list<std::string>& args)  {  return 1;  }
};

/********************************************************************************/

/** \brief Option that has one argument.
 *
 * This is a special option with only one argument.
 */
class OptionOneParam : public Option
{
public:

    /** Constructor.
     * \param[in] name : name of the option
     * \param[in] help : textual help for this option
     * \param[in] mandatory : tells whether this option is mandatory or not
     * \param[in] defaultValue : default value for the option
     * \param[in] multiple : tells whether this option may be used more than once
     * \param[in] include : list of names of options that must be used with the current one
     * \param[in] exclude : list of names of options that must not be used with the current one
     */
    OptionOneParam (
        const std::string& name,
        const std::string& help,
        bool mandatory = false,
        const std::string& defaultValue = "",
        int multiple = 0,
        const char* include = "",
        const char* exclude = ""
    )
        : Option (name, 1, mandatory, defaultValue, help, multiple, include, exclude)
    {
    }

    /** \copydoc Option::proceed */
    int proceed (const std::list<std::string>& args)
    {
        _param = *(args.begin());
        return 1;
    }
};

/********************************************************************************/


/** \brief Parser that analyzes command line options.
 *
 * Client can use this class for registering command line options specifications
 * and then can use it for parsing some command line options, typically given
 * as arguments of the 'main' function.
 *
 * Code sample:
 * \code
 * int main (int argc, char* argv[])
 * {
 *      // we create a parser
 *      OptionsParser parser;
 *
 *      // we register some options to it
 *      parser.add (new OptionOneParam ("-p", "Program Name [plastp, tplastn, plastx or tplastx]") );
 *      parser.add (new OptionOneParam ("-d", "Subject database file") );
 *      parser.add (new OptionOneParam ("-i", "Query database file") );
 *      parser.add (new OptionOneParam ("-h", "Help") );
 *
 *      // we parse the provided options
 *      int nbErrors = parser.parse (argc, argv);
 *
 *      // we retrieve options information as properties
 *      dp::IProperties* props = parser.getProperties ();
 * }
 * \endcode
 */
class OptionsParser : public system::SmartPointer
{
public:

    /** Constructor. */
    OptionsParser (const std::string& name="");

    /** Destructor. */
    ~OptionsParser ();

    /** Add the options of the provided parser.
     * \param[in] parser : the parser from which we add the options. */
    void push_back (const OptionsParser& parser);

    /** Add the options of the provided parser.
     * \param[in] parser : the parser from which we add the options. */
    void push_front (OptionsParser& parser);

    /** Add the options of the provided parser.
     * \param[in] parser : the parser from which we add the options. */
    void add (OptionsParser& parser) { push_back (parser); }

    /** Add an option to the pool of recognized options.
     * \param[in] opt : option to be registered to the parser.
     * \return the number of known options
     */
    int push_back (Option* opt);

	
	/** remove an option to the pool of recognized options.
     * \param[in] label : label to be removed to the parser.
     * \return the number of known options
     */
    int remove (const char * label);
	
	
    /** Add an option to the pool of recognized options.
     * \param[in] opt : option to be registered to the parser.
     * \return the number of known options
     */
    int push_front (Option* opt);

    /** Perform the analyze of the arguments.
     * \param[in] argc : number of command line arguments.
     * \param[in] argv : table of arguments
     * \return number of parsing errors.
     */
    misc::IProperties* parse (int argc, char* argv[]);

    /** Perform the analyze of the arguments.
     * \param[in] s : string containing the options to be parsed
     * \return number of parsing errors.
     */
    misc::IProperties*  parse (const std::string& s);

    /** Display errors (if there are some).
     * \param[in] fp : the file descriptor where to dump the errors
     */
    void displayErrors (FILE* fp = stdout);

    /** Display warnings (if there are some).
     * \param[in] fp : the file descriptor where to dump the warnings
     */
    void displayWarnings (FILE* fp = stdout);

    /** Display the help of each options recorded.
     * \param[in] fp : the file descriptor where to dump the help
     */
    void displayHelp (FILE* fp = stdout);

    /** Display version and other information
     * \param[in] fp : the file descriptor where to dump the information
     */
    void displayVersion (FILE* fp = stdout);

    /** Tells (after Proceed) if one option whose name is given has been seen or not.
     * \param[in] txt : the option name to be checked
     * \return true if option was seen, false otherwise.
     */
    bool saw (const std::string& txt);

    /** Return the list of seen options during the parsing.
     * \return the list of seen options.
     */
    const std::list<Option*>& getSeenOptions ()  { return _seenOptions; }

    /** Tells whether an option has been seen or not, given its label.
     * \return true if seed, false otherwise.
     */
    const Option* getSeenOption (const std::string& label);

    /** Return a IProperties instance holding parsed options information.
     * \return the IProperties instance.
     */
    misc::IProperties* getProperties ()  { return _properties; }

    /** Return the list of options that define the parser.
     * \return the list of options.
     */
    const std::list<Option*>& getOptions () const { return _options; }

    /** */
    void setName (const std::string& name)  { _name = name; }

private:

    /** */
    std::string _name;

    /** */
    misc::IProperties* _properties;
    void setProperties (misc::IProperties* properties)  { SP_SETATTR(properties); }
    void buildProperties ();

    /** List of Options*. */
    std::list<Option*> _options;

    /** List of errors. */
    std::list<std::string> _errors;

    /** List of Text* of warnings. */
    std::list<std::string> _warnings;

    /** List of seen options. */
    std::list<Option*> _seenOptions;

    /** */
    char _proceed;

    /** */
    int _argc;
    /** */
    char** _argv;
    /** */
    int _currentArg;

    /** */
    Option* lookForOption (char* txt);

    /** */
    char* nextArg ();

    /** */
    void getOptionArgs (const Option* option, std::list<std::string>& args);

    /** */
    void giveToNoOption (char* txt);

    /** */
    char* checkExcludingOptions (const Option* option);

    /** */
    void checkIncludingOptions ();

    /** */
    void checkMandatoryOptions ();
};

/********************************************************************************/

/** \brief Exception class to be used for option management error.
 *
 * This class should be thrown when something went wrong during options parsing.
 */
class OptionFailure
{
public:

    /** Constructor.
     * \param[in] parser : the parser that threw the exception.
     */
    OptionFailure (OptionsParser& parser) : _parser(parser)  {}

    /** Getter on the parser.
     * \return the parser.
     */
    OptionsParser& getParser ()  { return _parser; }

private:
    OptionsParser& _parser;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_OPTION_PARSER_HPP_ */
