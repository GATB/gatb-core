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

#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <memory>

using namespace std;

#define ABS(a)  ((a)<0 ? -(a) : (a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/********************************************************************************/

/** \brief Test class for operating system operations
 */
class TestIterable : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestIterable);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestIterable);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestIterable);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
