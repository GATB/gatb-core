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

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_ALGORITHM_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
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

/** \brief Framework class for implementing algorithm
 */
class Algorithm : public dp::SmartPointer
{
public:

    /** Constructor.
     * \param[in] name: name of the algorithm. */
    Algorithm (const std::string& name, gatb::core::tools::misc::IProperties* input=0);

    /** Destructor. */
    virtual ~Algorithm ();

    /** Get tool name
     * \return the algorithm name. */
    std::string getName () const  { return _name; }

    /** */
    virtual void execute () = 0;

    /** */
    virtual IProperties*            getInput      ()  { return _input;      }
    virtual IProperties*            getOutput     ()  { return _output;     }
    virtual IProperties*            getInfo       ()  { return _info;       }
    virtual dp::ICommandDispatcher* getDispatcher ()  { return _dispatcher; }
    virtual TimeInfo&               getTimeInfo   ()  { return _timeInfo;   }

    /** */
    template<typename Item> dp::Iterator<Item>* createIterator (dp::Iterator<Item>* iter, size_t nbIterations=0, const char* message=0)
    {
        if (nbIterations > 0 && message != 0)
        {
            //  We create some listener to be notified every 1000 iterations and attach it to the iterator.
            dp::impl::SubjectIterator<Item>* iterSubject = new dp::impl::SubjectIterator<Item> (iter, nbIterations/100);
            iterSubject->addObserver (createIteratorListener (nbIterations, message));

            /** We assign the used iterator to be the subject iterator. */
            iter = iterSubject;
        }

        /** We return the result. */
        return iter;
    }

    /** */
    virtual dp::IteratorListener* createIteratorListener (size_t nbIterations, const char* message);

protected:

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUriByKey (const std::string& key)  { return getUri (getInput()->getStr(key)); }

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUri (const std::string& str)  { return getInput()->getStr(STR_PREFIX) + str; }

    /** Setters. */
    void setInput      (IProperties*            input)       { SP_SETATTR (input);      }
    void setOutput     (IProperties*            output)      { SP_SETATTR (output);     }
    void setInfo       (IProperties*            info)        { SP_SETATTR (info);       }
    void setDispatcher (dp::ICommandDispatcher* dispatcher)  { SP_SETATTR (dispatcher); }

private:

    /** Name of the tool (set at construction). */
    std::string _name;

    IProperties* _input;

    IProperties* _output;

    IProperties* _info;

    dp::ICommandDispatcher* _dispatcher;

    /** */
    TimeInfo _timeInfo;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_ALGORITHM_HPP_ */
