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

#include <CppunitCommon.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/system/impl/TimeCommon.hpp>

#include <gatb/tools/misc/api/Range.hpp>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <memory>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::misc;

#define ABS(a)  ((a)<0 ? -(a) : (a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for miscellaneous operations
 */
class TestMisc : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestMisc);

        CPPUNIT_TEST_GATB (range_checkIterator1);
        CPPUNIT_TEST_GATB (range_checkIterator2);
        //  NEED A GOOD TIMER ACCURACY...  CPPUNIT_TEST_GATB (range_checkPerformance);

    CPPUNIT_TEST_SUITE_GATB_END();

public:
    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/
    /** \brief test of Range class and its associated iterator.
     *
     * We create an iterator for an integer range and check that each iterated item is ok.
     */
    void range_checkIterator1 ()
    {
        /** We create a range. */
        Range<size_t> range (1, 1000);

        /** We create an iterator from the range. */
        Range<size_t>::Iterator it (range);

        size_t check = range.getBegin();

        /** We iterate each item of the range. */
        for (it.first(); !it.isDone(); it.next())
        {
            /** We check that the current iterated item is ok. */
            CPPUNIT_ASSERT (it.item() == check++);
        }
    }

    /********************************************************************************/
    /** \brief test of Range associated iterator.
     *
     * We create an iterator for an integer range and check that each iterated item is ok.
     */
    void range_checkIterator2 ()
    {
        size_t from=1, to=100, check=from;

        /** We create an iterator from the range. */
        Range<size_t>::Iterator it (from, to);

        /** We iterate each item of the range. */
        for (it.first(); !it.isDone(); it.next())
        {
            /** We check that the current iterated item is ok. */
            CPPUNIT_ASSERT (it.item() == check++);
        }
    }

    /********************************************************************************/
    void range_checkPerformance ()
    {
        typedef u_int64_t INT;

        u_int64_t sum1=0, sum2=0;

        INT begin=1, end=1000*1000*1000;

        /** We create a range. */
        Range<INT> range (begin, end);

        /** We create an iterator from the range (seen as an Iterable). */
        Range<INT>::Iterator it (range);

        /** We need something to measure elapsed time. */
        ITime& ts = System::time();

        /** Note that we use 'volatile' keyword in order to avoid unwanted optimizations here. */
        volatile ITime::Value t0 = ts.getTimeStamp();

        for (it.first(); !it.isDone(); it.next())   {  sum1 += it.item();  }

        volatile ITime::Value t1 = ts.getTimeStamp();

        for (INT i=begin; i<=end; i++)   {  sum2 += i;  }

        volatile ITime::Value t2 = ts.getTimeStamp();

        /** We check we got the same result with both loops. */
        CPPUNIT_ASSERT (sum1 > 0  &&  sum1 == sum2);

        /** We check that performances are of same order. */
        double err = 100.0 * ABS (1.0 - (double)(t1-t0) / (double)(t2-t1));

        /** The 'iterator' loop should be less than 100% slower than the 'direct' loop. */
        CPPUNIT_ASSERT (err < 100);
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION (TestMisc);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
