/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Exception.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Operating System common abstraction.
 */

#ifndef _GATB_CORE_SYSTEM_IRESOURCE_HPP_
#define _GATB_CORE_SYSTEM_IRESOURCE_HPP_

/********************************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <string>

/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writing
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Operating System abstraction layer */
namespace system    {
/********************************************************************************/

 /** \brief Exception class for operating system failures
  */
 class Exception
 {
 public:

     /** Default constructor. */
     Exception ()  {}

     /** Constructor with some information in a "printf" way.
      * \param[in] format : format of the message
      * \param[in] ... : arguments for the format parameter
      */
     Exception (const char* format, ...)
     {
         va_list args;  va_start (args, format);  init (format, args);  va_end (args);
     }

     /** Returns the description message.
      * \return the message
      */
     const char* getMessage ()  { return _message.c_str(); }

 protected:

     /** */
     void init (const char* format, va_list args)
     {
         char buffer[256];
         vsnprintf (buffer, sizeof(buffer), format, args);
         _message.assign (buffer);
     }

     /** The informative message. */
     std::string _message;
 };

 /********************************************************************************/

 /** \brief Exception class with information got from strerror_r
  */
 class ExceptionErrno : public Exception
 {
 public:

    /** Constructor. The error message is built by calling strerror function. */
    ExceptionErrno (const char* format, ...)
    {
        va_list args;  va_start (args, format);  init (format, args);  va_end (args);

        char* buffer = (char*) malloc (BUFSIZ);
        if (buffer != NULL)
        {
            *buffer = 0;

            strerror_r (errno, buffer, BUFSIZ);
            {  _message += std::string(" (") + std::string(buffer) + std::string(")");  }
            free(buffer);
        }
    }
 };

 /********************************************************************************/

 /** \brief Exception class for lack of implementation
  */
 class ExceptionNotImplemented : public Exception
 {
 public:

     /** Constructor. */
     ExceptionNotImplemented ()  {  _message = "NOT IMPLEMENTED";  }
 };

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IRESOURCE_HPP_ */
