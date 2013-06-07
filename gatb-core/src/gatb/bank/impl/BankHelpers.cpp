/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <iostream>

using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperties* BankHelper::convert (IBank& in, IBank& out, tools::dp::impl::IteratorListener* progress)
{
    // We need an iterator on the input bank.
    Iterator<Sequence>* itBank = in.iterator();
    LOCAL (itBank);

    SubjectIterator<Sequence> itSeq (itBank, 100*1000);

    if (progress != 0)  {  itSeq.addObserver (progress);  }

    u_int64_t   nbSeq = 0;
    u_int64_t sizeSeq = 0;

    // We get current time stamp
    TimeSystem ts (ITime::MSEC);
    ITime::Value t0 = ts.getTimeStamp();

    for (itSeq.first(); !itSeq.isDone(); itSeq.next())
    {
        nbSeq ++;
        sizeSeq += (itSeq)->getDataSize();

        out.insert ( itSeq.item());
    }

    /** We flush the output bank (important if it is a file output). */
    out.flush ();

    // We get current time stamp
    ITime::Value t1 = ts.getTimeStamp();

    /** We create the properties result. */
    Properties* props = new Properties ();

    props->add (0, "conversion", "");

    props->add (1, "time_sec",         "%.2f", (double)(t1-t0)/1000.0);
    props->add (1, "sequences_number", "%ld",  nbSeq);
    props->add (1, "sequences_size",   "%ld",  sizeSeq);
    props->add (1, "output_size",      "%ld",  out.getSize());

    /** We return the result.*/
    return props;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
