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

#include <iostream>

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
    return dbprefix +  a;
}

/********************************************************************************/
int main (int argc, char **argv)
{
    /** We set the directory where the db are. */
    dbprefix    = (argc >=2 ? argv[1] : "../test/db/");

    /** We may have an ouput xml file. */
    char* xmloutput = (argc >=3 ? argv[2] : 0);

    // informs test-listener about testresults
    TestResult testresult;

    // register listener for collecting the test-results
    TestResultCollector collectedresults;
    testresult.addListener (&collectedresults);

#if 0
    BriefTestProgressListener progress;
    testresult.addListener (&progress);
#endif

    TextTestRunner runner;
    runner.addTest ( TestFactoryRegistry::getRegistry().makeTest ());
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

    // return 0 if tests were successful
    return collectedresults.wasSuccessful() ? 0 : 1;
}
