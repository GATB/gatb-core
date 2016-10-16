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

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/collections/impl/CollectionCache.hpp>

#include <gatb/tools/misc/api/Range.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>

#define USE_LARGEINT_CONSTRUCTOR 1 // got lazy to change the unit test to avoid construction.
#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/LargeInt.hpp>

using namespace std;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;

using namespace gatb::core::tools::math;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

#define ABS(a)  ((a)<0 ? -(a) : (a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for operating system operations
 */
class TestCollection : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestCollection);

        CPPUNIT_TEST_GATB (collection_check1);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    template<typename T> void collection_check1_aux (T* table, size_t nb)
    {
        // new since STORAGE_FILE modification in october 2016: remove the file, because BagFile doesn't anymore
         System::file().remove ("foo");

        /** We create a collection. */
        Collection<T>* c = new CollectionFile<T> ("foo");
        LOCAL(c);

        /** We insert some items. */
        for (size_t i=0; i<nb; i++)  { c->insert (table[i]); }

        /** We have to flush to be sure all items are in the collection. */
        c->flush();

        if (c->getNbItems() != (int)nb)
        {
            std::cout << "in anticipation of unit test failing, c->getNbItems() = " << c->getNbItems() << " , nb = " << nb << std::endl;
        }

        /** We check we have the correct number of items. */
        CPPUNIT_ASSERT (c->getNbItems() == (int)nb);

        /** We iterate the collection. */
        Iterator<T>* it = c->iterator();
        LOCAL(it);

        /** We check we get the correct items. */
        size_t i=0;
        for (it->first(); !it->isDone(); it->next())  {  CPPUNIT_ASSERT (it->item() == table[i++]);  }

        /** We delete physically the collection. */
        //c->remove ();
    }

    /********************************************************************************/
    void collection_check1 ()
    {
        int table1[] = { 1,2,3,5,8,13,21,34};
        collection_check1_aux<int> (table1, ARRAY_SIZE(table1));

        string table2[] = { "a", "bc", "def", "ghij"};
        collection_check1_aux<string> (table2, ARRAY_SIZE(table2));

        // Test table3[] = { Test(1, 1025, 3.14), Test(3, 298, 2.713), Test(19, 9874, 0.577)};
        // collection_check1_aux<Test> (table3, ARRAY_SIZE(table3));

        NativeInt64 table4[] = { 12354684684, 6876436549, 87654351, 6843516877, 68435434874};
        collection_check1_aux<NativeInt64> (table4, ARRAY_SIZE(table4));

        LargeInt<3> table5[] = { LargeInt<3>(413434), LargeInt<3>(987654123), LargeInt<3>(123), LargeInt<3>(1) };
        collection_check1_aux<LargeInt<3> > (table5, ARRAY_SIZE(table5));
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestCollection);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestCollection);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
