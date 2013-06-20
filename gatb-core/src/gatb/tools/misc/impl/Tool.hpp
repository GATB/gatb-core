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
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>

#include <string>
#include <list>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Framework class for implementing tools (ie. binary tools).
 */
class Tool : public dp::SmartPointer
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
    virtual IProperties* run (IProperties* input);

    /** */
    virtual IProperties* run (int argc, char* argv[]);

    /** */
    virtual IProperties*            getInput      ()  { return _input;      }
    virtual IProperties*            getOutput     ()  { return _output;     }
    virtual IProperties*            getInfo       ()  { return _info;       }
    virtual OptionsParser*          getParser     ()  { return _parser;     }
    virtual dp::ICommandDispatcher* getDispatcher ()  { return _dispatcher; }
    virtual TimeInfo&               getTimeInfo   ()  { return _timeInfo;   }

    /** */
    template<typename Item> dp::Iterator<Item>* createIterator (dp::Iterator<Item>* iter, size_t nbIterations=0, const char* message=0)
    {
        if (_input->get(STR_QUIET) == 0 && _input->get(STR_NO_PROG_BAR) == 0)
        {
            if (nbIterations > 0 && message != 0)
            {
                //  We create some listener to be notified every 1000 iterations and attach it to the iterator.
                dp::impl::SubjectIterator<Item>* iterSubject = new dp::impl::SubjectIterator<Item> (iter, nbIterations/100);
                iterSubject->addObserver (createIteratorListener (nbIterations, message));

                /** We assign the used iterator to be the subject iterator. */
                iter = iterSubject;
            }
        }

        /** We return the result. */
        return iter;
    }

    /** */
    static const char* STR_NB_CORES;
    static const char* STR_STATS_XML;
    static const char* STR_QUIET;
    static const char* STR_NO_PROG_BAR;
    static const char* STR_URI_PREFIX;
    static const char* STR_URI_DATABASE;
    static const char* STR_URI_OUTPUT;
    static const char* STR_HELP;

protected:

    /** */
    virtual void execute () = 0;

    /** */
    virtual void preExecute  ();
    virtual void postExecute ();

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUriByKey (const std::string& key)  { return getUri (_input->getStr(key)); }

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUri (const std::string& str)  { return _input->getStr(STR_URI_PREFIX) + str; }

    /** Setters. */
    void setInput      (IProperties*            input)       { SP_SETATTR (input);      }
    void setOutput     (IProperties*            output)      { SP_SETATTR (output);     }
    void setInfo       (IProperties*            info)        { SP_SETATTR (info);       }
    void setParser     (OptionsParser*          parser)      { SP_SETATTR (parser);     }
    void setDispatcher (dp::ICommandDispatcher* dispatcher)  { SP_SETATTR (dispatcher); }

private:

    /** Name of the tool (set at construction). */
    std::string _name;

    IProperties* _input;

    IProperties* _output;

    IProperties* _info;

    OptionsParser* _parser;

    dp::ICommandDispatcher* _dispatcher;

    /** */
    TimeInfo _timeInfo;

    /** */
    virtual dp::IteratorListener* createIteratorListener (size_t nbIterations, const char* message);

    friend class ToolComposite;
};

/********************************************************************************/

class ToolComposite : public Tool
{
public:

    /** Constructor.
     * \param[in] name: name of the tool. */
    ToolComposite (const std::string& name = "tool");

    /** */
    ~ToolComposite ();

    /** */
    IProperties* run (int argc, char* argv[]);

    /** */
    void add (Tool* tool);

private:

    std::list<Tool*> _tools;

    /** */
    void execute ();
    void preExecute  ();
    void postExecute ();
};

/********************************************************************************/

/** */
class ToolProxy : public Tool
{
public:

    /** */
    ToolProxy (Tool* ref) : Tool("proxy"), _ref (ref)  {}

    /** */
    virtual OptionsParser* getParser ()  {  return _ref->getParser();  }

    /** */
    virtual IProperties* getInput  ()  { return _ref->getInput();     }
    virtual IProperties* getOutput ()  { return _ref->getOutput();    }
    virtual IProperties* getInfo   ()  { return _ref->getInfo();      }

    /** */
    virtual dp::ICommandDispatcher*    getDispatcher ()   { return _ref->getDispatcher(); }

    /** */
    virtual TimeInfo&    getTimeInfo ()   { return _ref->getTimeInfo (); }

    /** */
    Tool* getRef ()  { return _ref; }

private:
    Tool* _ref;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_ */
