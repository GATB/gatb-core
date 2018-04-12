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

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/api/Sequence.hpp>


#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;


extern std::string DBPATH (const string& a);

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for the Sequence object
 */
class TestSequence : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestSequence);

    CPPUNIT_TEST_GATB(bank_checkRevcomp);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/

    void bank_checkRevcomp ()
    {
        char* buf = (char*)"ACTACGATCGATGTA";
        Sequence s1 (buf);

        string t;
        t = s1.getRevcomp();
		CPPUNIT_ASSERT(t.compare("TACATCGATCGTAGT")==0);



        char* buf2 = (char*)"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA";
        Sequence s2 (buf2);

        t = s2.getRevcomp();
		CPPUNIT_ASSERT(t.compare("TTGGTTAAAGTATTTAGTGACCTAAGTCAATAAAATTTTAATTTACTCACGGCAGGTAACCAGTTCAGAA")==0);

    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestSequence);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestSequence);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

