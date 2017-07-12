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
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/collections/impl/MapMPHF.hpp>
#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>
#include <gatb/system/api/Exception.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <memory>

using namespace std;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;
using namespace gatb::core::tools::math;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

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
        CPPUNIT_TEST_GATB (checkMapMPHF);

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

        T badKey;
        badKey.setVal(hash.getMaxNbItems() + 100);

        /** We insert the maximum number of items. */
        for (int i=1; i<=hash.getMaxNbItems(); i++)  {  
            T idx; idx.setVal(i);
            CPPUNIT_ASSERT_NO_THROW (hash.increment (idx));  
        }

        /** We add a new key => we should get an exception. */
        CPPUNIT_ASSERT_THROW (hash.increment (badKey), core::system::Exception);

        /** We check that we have all the required keys. */
        for (int i=1; i<=hash.getMaxNbItems(); i++)  {  
            T idx; idx.setVal(i);
            CPPUNIT_ASSERT (hash.get (idx) == true);  
        }

        /** We check that we don't have an non registered key. */
        CPPUNIT_ASSERT (hash.get (badKey) == false);

        /** We iterate the map. */
        Iterator <Abundance<T> >* it = hash.iterator();
        LOCAL (it);

        size_t nbItems = 0;
        for (it->first(); !it->isDone(); it->next(), nbItems++)
        {
            /** All abundances should be one. */
            CPPUNIT_ASSERT (it->item().abundance == 1);
        }

        CPPUNIT_ASSERT ((int)nbItems == hash.getMaxNbItems());
    }

    /********************************************************************************/
    void checkOAHash ()
    {
        size_t table[] = { 1024, 10*1024, 100*1024, 1000*1024};

        for (size_t i=0; i<ARRAY_SIZE(table); i++)
        {
            checkOAHash_aux<NativeInt64>  (table[i]);
    #if INT128_FOUND == 1
            checkOAHash_aux<NativeInt128> (table[i]);
    #endif
            checkOAHash_aux<LargeInt<3> >  (table[i]);
            checkOAHash_aux<LargeInt<4> >  (table[i]);
            checkOAHash_aux<LargeInt<5> >  (table[i]);
        }
    }

    /********************************************************************************/
    static void checkMapMPHF_progress (size_t round, size_t initial, size_t remaining)
    {
    }

    /** */
    void checkMapMPHF ()
    {
        float val = 0;

        u_int8_t keysValue[] = {14, 35, 1, 9, 65, 37, 12, 24, 98, 124, 32};

        /** We create a file with some keys (coded as NativeInt8 values). */
        const char* filename = "keys";
        System::file().remove (filename); // new since october 2016 BagFile change
        BagFile<NativeInt8> keysFile (filename);
        for (size_t i=0; i<ARRAY_SIZE(keysValue); i++)  {  keysFile.insert (keysValue[i]); }
        keysFile.flush();

        /** We create some Iterable for our keys. */
        IterableFile<NativeInt8> keys(filename);

        CPPUNIT_ASSERT (keys.getNbItems() == ARRAY_SIZE(keysValue));

        /** We create a map for our keys; the values are floats. */
        MapMPHF <NativeInt8, float> map1;

        /** We build the hash function for the given keys. */
        map1.build (keys);

        CPPUNIT_ASSERT (map1.size() == (size_t)keys.getNbItems());

        /** We populate the values. */
        Iterator<NativeInt8>* itKeys = keys.iterator();   LOCAL (itKeys);
        for (itKeys->first(); !itKeys->isDone(); itKeys->next(), val++)
        {
            /** We change the value for the current key. */
            map1.at(itKeys->item()) = val;
        }

        /** We check the values. */
        val = 0;
        for (itKeys->first(); !itKeys->isDone(); itKeys->next(), val++)
        {
            /** We check the value for the current key. */
            CPPUNIT_ASSERT (map1.at(itKeys->item()) == val);
        }

        /** We create a storage object. */
        Storage* storage = StorageFactory(STORAGE_HDF5).create ("map", true, false);
        LOCAL (storage);

        /** We save the map. */
        map1.save (storage->root(), "mphf");

        /** We create another map. */
        MapMPHF <NativeInt8, float> map2;

        /** We load the hash function. */
        map2.load (storage->root(), "mphf");

        CPPUNIT_ASSERT (map2.size() == (size_t)keys.getNbItems());

        /** We populate the values. */
        val = 0;
        for (itKeys->first(); !itKeys->isDone(); itKeys->next(), val++)
        {
            /** We change the value for the current key. */
            map2.at(itKeys->item()) = val;
        }

        /** We check the values. */
        val = 0;
        for (itKeys->first(); !itKeys->isDone(); itKeys->next(), val++)
        {
            /** We check the value for the current key. */
            CPPUNIT_ASSERT (map2.at(itKeys->item()) == val);
        }

        /** We compare the values of the two maps. */
        for (itKeys->first(); !itKeys->isDone(); itKeys->next(), val++)
        {
            /** We check the value for the current key. */
            CPPUNIT_ASSERT (map1.at(itKeys->item()) == map2.at(itKeys->item()));
        }

        /** We compare the values of the two maps (index iteration) */
        CPPUNIT_ASSERT (map1.size() == map2.size());
        for (size_t i=0; i<map1.size(); i++)
        {
            /** We check the value for the current key. */
            CPPUNIT_ASSERT (map1.at(i) == map2.at(i));
        }
    }
};

/********************************************************************************/

//Temporary desactivate TestMap (for continuous integration purpose)
CPPUNIT_TEST_SUITE_REGISTRATION      (TestMap);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestMap);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
