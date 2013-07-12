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

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/system/api/Exception.hpp>

#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <memory>

using namespace std;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;
using namespace gatb::core::tools::math;
using namespace gatb::core::system;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/********************************************************************************/

/** \brief Test class for operating system operations
 */
class TestMap : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestMap);

        CPPUNIT_TEST_GATB (checkOAHash);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/
    template<typename T>
    void checkOAHash_aux (size_t maxMemory)
    {
        /** We create a hash with a maximum memory size. */
        OAHash <T> hash (maxMemory);

        size_t badKey = hash.getMaxNbItems() + 100;

        /** We insert the maximum number of items. */
        for (size_t i=1; i<=hash.getMaxNbItems(); i++)  {  CPPUNIT_ASSERT_NO_THROW (hash.increment (i));  }

        /** We add a new key => we should get an exception. */
        CPPUNIT_ASSERT_THROW (hash.increment (badKey), core::system::Exception);

        /** We check that we have all the required keys. */
        for (size_t i=1; i<=hash.getMaxNbItems(); i++)  {  CPPUNIT_ASSERT (hash.get (i) == true);  }

        /** We check that we don't have an non registered key. */
        CPPUNIT_ASSERT (hash.get (badKey) == false);

        /** We iterate the map. */
        Iterator <pair<T,u_int32_t> >* it = hash.iterator();
        LOCAL (it);

        size_t nbItems = 0;
        for (it->first(); !it->isDone(); it->next(), nbItems++)
        {
            /** All abundances should be one. */
            CPPUNIT_ASSERT (it->item().second == 1);
        }

        CPPUNIT_ASSERT (nbItems == hash.getMaxNbItems());
    }

    /********************************************************************************/
    void checkOAHash ()
    {
        size_t table[] = { 1024, 10*1024, 100*1024, 1000*1024};

        for (size_t i=0; i<ARRAY_SIZE(table); i++)
        {
            checkOAHash_aux<NativeInt64>  (table[i]);
    #ifdef INT128_FOUND
            checkOAHash_aux<NativeInt128> (table[i]);
    #endif
            checkOAHash_aux<LargeInt<3> >  (table[i]);
            checkOAHash_aux<LargeInt<4> >  (table[i]);
            checkOAHash_aux<LargeInt<5> >  (table[i]);
        }
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION (TestMap);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
