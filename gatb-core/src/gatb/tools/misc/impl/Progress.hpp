/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
 *****************************************************************************/

/** \file Progress.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Progress feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_PROGRESS_HPP_
#define _GATB_CORE_TOOLS_MISC_PROGRESS_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <string>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Progress information display feature.
 *
 * This class is used to display some progression information, for instance during
 * the iteration of an iterator.
 *
 * It provides 3 methods that are supposed to be called by the progressing job:
 *      - init   : called at the beginning of the job
 *      - finish : called at the end of the job
 *      - inc    : called to notify some amount of job being done
 *
 * It can be used as this by the SubjectIterator class that encapsulates an iterator
 * and calls some functor for notifying progress information.
 *
 * This implementation merely dumps '-' characters as progression goes on.
 *
 * \see gatb::tools::dp::impl::SubjectIterator
 */
class Progress : public dp::impl::IteratorListener
{
public:

    /** Constructor.
     * \param[in] ntasks : nb of items to be processed
     * \param[in] msg : message to be displayed
     * \param[in] os  : output stream where progress information is to be displayed
     */
    Progress (u_int64_t ntasks, const char * msg, std::ostream& os=std::cerr);

    /** Destructor. */
    virtual ~Progress () {}

    /** Initialization of the object. */
    void init ();

    /** Finish the progress information. */
    void finish ();

    /** Increase the number of currently done tasks.
     * \param[in] ntasks_done : amount of job done before previous call. */
    void inc (u_int64_t ntasks_done);

    /** Set the current number of tasks done.
     * \param[in] ntasks_done :  sets the current number of job done. */
    void set (u_int64_t ntasks_done);

    /** \copydoc dp::impl::IteratorListener::setMessage*/
    void setMessage (const char* format, ...);

protected:

    virtual void update     ();
    virtual void postInit   ();
    virtual void postFinish ();

    std::string   message;
    u_int64_t     done;
    u_int64_t     todo;
    int           subdiv ; // progress printed every 1/subdiv of total to do
    double        partial;
    double        steps ; //steps = _todo/subidv
    std::ostream& os;
    char          buffer[256];
};

/********************************************************************************/

/** \brief Progress information display feature with timing information
 *
 * This is another way of dumping information progressing, that provides information
 * such as the elapsed time and an estimation of the remaining time.
 */
class  ProgressTimer : public Progress
{
public:

    /** Constructor.
     * \param[in] ntasks : nb of items to be processed
     * \param[in] msg : message to be displayed
     * \param[in] os  : output stream where progress information is to be displayed
     */
    ProgressTimer (u_int64_t ntasks, const char* msg, std::ostream& os=std::cerr);

protected:

    void update ();
    void postInit ();
    void postFinish ();

    system::ITime::Value  heure_debut;
    system::ITime::Value  heure_actuelle;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_PROGRESS_HPP_ */
