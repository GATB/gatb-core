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

#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/bcalm2/unionFind.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <iostream>
#include <thread>

using namespace std;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::collections::impl;

extern std::string DBPATH (const string& a);

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for genomic databases management
 */
class TestBcalm : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestBcalm);

        CPPUNIT_TEST_GATB (bcalm_test1); 
        CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}


    /********************************************************************************/
    void bcalm_test1 () // i wanna test de UF because there was a weird behavior in bcalm. turns out: it wasn't the UF that was faulty, but let's keep this test.
    {
        int nb_uf_elts = 3000000;
        unionFind uf(nb_uf_elts);


        auto doJoins = [&uf](int start, int end)
        {
            for (int i = start; i < end; i++)
            {
                uf.union_(i,i+1);
                //std::cout << "joining " << i << " " << i+1 << std::endl; // it's indeed threaded, this checks it
            }

        };

        doJoins(nb_uf_elts/2,nb_uf_elts-1);
        doJoins(0,nb_uf_elts/2);
        
        // what if it's not a multithread bug but an order of operations bug, somehow
        /* // threaded
        std::thread first(doJoins,0,nb_uf_elts/3);
        std::thread second(doJoins,nb_uf_elts/3,2*(nb_uf_elts/3));
        std::thread third(doJoins,2*(nb_uf_elts/3),nb_uf_elts-1);

        first.join();
        second.join();
        third.join();
        */

        int foundclass = uf.find(0); // not always 0, depends which thread starts first
        for (int i = 1; i < nb_uf_elts; i++)
        {   
            int ufclass = uf.find(i);
            if (ufclass != foundclass)
                std::cout<< "\nclass of " << i << " is " << ufclass <<std::endl;
            CPPUNIT_ASSERT (ufclass == foundclass);
        }

    }

};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestBcalm);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestBcalm);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

