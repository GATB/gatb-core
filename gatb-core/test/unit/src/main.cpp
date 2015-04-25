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

#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;

/********************************************************************************/
/** GATB package */
namespace gatb    {
/** GATB tests package */
namespace tests   {
} }
/********************************************************************************/

static string dbprefix;

/********************************************************************************/
std::string DBPATH (const string& a)
{
    return dbprefix + string("/") +  a;
}

/********************************************************************************/
int main (int argc, char **argv)
{
    /** We may launch only one test. */
    char* testname = strdup (argc >=2 ? argv[1] : "All Tests");

    /** We set the directory where the db are. */
    dbprefix = (argc >=3 ? argv[2] : "../test/db");

    /** We may have an ouput xml file. */
    char* xmloutput = (argc >=4 ? argv[3] : 0);

    // informs test-listener about testresults
    TestResult testresult;

    // register listener for collecting the test-results
    TestResultCollector collectedresults;
    testresult.addListener (&collectedresults);

#if 1
    BriefTestProgressListener progress;
    if (getenv ("CPPUNIT_VERBOSE")) {  testresult.addListener (&progress);  }
#endif

    TextTestRunner runner;

    for (char* loop = strtok (testname,","); loop != 0;  loop = strtok (NULL, ","))
    {
        runner.addTest ( TestFactoryRegistry::getRegistry(loop).makeTest ());
    }

    runner.run (testresult);

    // output results in compiler-format
    CompilerOutputter compileroutputter (&collectedresults, std::cout);
    compileroutputter.write ();

    // Output XML
    if (xmloutput != 0)
    {
        ofstream xmlFileOut (xmloutput);
        XmlOutputter xmlOut (&collectedresults, xmlFileOut);
        xmlOut.write();
    }

    free (testname);

    // return 0 if tests were successful
    return collectedresults.wasSuccessful() ? 0 : 1;
}
