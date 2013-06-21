/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/tools/misc/impl/Histogram.hpp>

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
void Histogram::inc (size_t abundance)
{
    _histogram [(abundance >= _length) ? _length : abundance] ++;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Histogram::save ()
{
    /** We first delete the previous file. */
    System::file().remove (_uri);

    IFile* file = System::file().newFile (_uri, "w");
    if (file != 0)
    {
        for (size_t cc=1; cc<_histogram.size(); cc++)  {  file->print ("%i\t%u\n", cc, _histogram[cc]);  }
        delete file;
    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
