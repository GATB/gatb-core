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

#ifndef  _GATB_CPPUNIT_COMMON_H_
#define  _GATB_CPPUNIT_COMMON_H_

/********************************************************************************/

#include <cppunit/extensions/TestFactoryRegistry.h>

#include <cppunit/TestAssert.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>

#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>

#include <cppunit/ui/text/TestRunner.h>

#include <cppunit/extensions/HelperMacros.h>

/********************************************************************************/

/* NOTE ABOUT THE MACROS ABOVE...
 * If we use default CPPUNIT_TEST_SUITE and CPPUNIT_TEST macros, having our test classes in namespaces,
 * we will have XML output files with mangled names, which is not fine.
 *
 * Note that we still want to have the test classes because we generate documentation from tests and
 * we hope that they will be in a specific location with all the library documentation.
 *
 * One solution is to define specific macros for GATB that avoids to use the RTTI system used by CppUnit,
 * which may lead in some implementation to mangled C++ names.
 */

#define CPPUNIT_TEST_SUITE_GATB(clazz)                              \
        static std::string getClazz() { return std::string(#clazz); }  \
        CPPUNIT_TEST_SUITE (clazz)

#define CPPUNIT_TEST_GATB(testMethod )                \
    CPPUNIT_TEST_SUITE_ADD_TEST(                            \
        ( new CPPUNIT_NS::TestCaller<TestFixtureType>(      \
                getClazz() + std::string("::") + std::string(#testMethod),  \
                  &TestFixtureType::testMethod,             \
                  context.makeFixture() ) ) )

#define CPPUNIT_TEST_SUITE_GATB_END  CPPUNIT_TEST_SUITE_END


#define CPPUNIT_TEST_SUITE_REGISTRATION_GATB(suite)         \
    CPPUNIT_TEST_SUITE_NAMED_REGISTRATION (suite,#suite)

using namespace CppUnit;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Abstract class for unit testing.
 *
 * Just here for having our own factorizing name (instead of CPPUNIT TestFixture name).
 */
class Test : public TestFixture
{
};

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

/********************************************************************************/


#endif /* _GATB_CPPUNIT_COMMON_H_ */
