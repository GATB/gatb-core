/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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
