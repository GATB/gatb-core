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

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/api/Macros.hpp>
#include <gatb/tools/misc/api/Range.hpp>

#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/system/impl/System.hpp>

#include <vector>

using namespace std;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

#define ABS(a)  ((a)<0 ? -(a) : (a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for operating system operations
 */
class TestBag : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestBag);

        CPPUNIT_TEST_GATB (bag_checkFile);
        CPPUNIT_TEST_GATB (bag_checkCache);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    void bag_checkFile ()
    {
        u_int32_t values[] = { 1, 2, 3, 5, 8, 13, 21, 34, 55, 89};
        size_t nbItems = ARRAY_SIZE(values);

        /** We build the uri of the test file. */
        string filename = System::file().getTemporaryDirectory() + string("/bagfile");

        /** We instantiate a bag. */
        BagFile<u_int32_t> bag (filename);

        /** We write some values into the bag. */
        for (u_int32_t i=0; i<nbItems; i++)  {  bag.insert (values[i]);  }

        /** We flush the bag to be sure that all items are put into the file. */
        bag.flush();

        /** We check the size of the file. */
        CPPUNIT_ASSERT (System::file().getSize(filename) == nbItems * sizeof (u_int32_t));

        /** We read the file. */
        size_t nbItemsRead = 0;
        IteratorFile<u_int32_t> it (filename);
        for (it.first(); !it.isDone(); it.next(), nbItemsRead++)
        {
            /** We check the current read item. */
            CPPUNIT_ASSERT (*it == values[nbItemsRead]);
        }

        /** We check we read the correct number of items. */
        CPPUNIT_ASSERT (nbItemsRead == nbItems);

        /** We remove the temporary file. */
        System::file().remove (filename);
    }

    /********************************************************************************/

    template <typename Item> class BagInsertCommand : public ICommand, public SmartPointer
    {
    public:
        BagInsertCommand (Bag<Item>* ref, size_t cacheSize, ISynchronizer* synchro, size_t nbIters)
            : _ref(ref), _cache(_ref, cacheSize, synchro), _nbIters(nbIters) {}

        void execute ()
        {
            for (u_int32_t i=1; i<=_nbIters; i++)  {  _cache.insert (i);  }

            _cache.flush();
        }

    private:
        Bag<Item>*     _ref;
        BagCache<Item> _cache;
        size_t         _nbIters;
    };

    /** */
    void bag_checkCache_aux (size_t nbItersTotal, size_t nbCores, size_t cacheSize)
    {
        /** We build the uri of the test file. */
        string filename = System::file().getTemporaryDirectory() + string("/bagfile");

        /** We set the number of iterations to be done in each thread. */
        u_int64_t nbItersPerThread = nbItersTotal / nbCores;

        /** We instantiate a bag. */
        BagFile<u_int32_t>* bag = new BagFile<u_int32_t>(filename);
        LOCAL (bag);

        /** We need a command dispatcher. */
        Dispatcher dispatcher (nbCores);

        /** We create several commands that will insert items into the bag file through a cache. */
        vector<ICommand*> commands (dispatcher.getExecutionUnitsNumber());
        for (size_t i=0; i<commands.size(); i++)  { commands[i] = new BagInsertCommand<u_int32_t> (bag, cacheSize, System::thread().newSynchronizer(), nbItersPerThread); }

        /** We dispatch the commands. */
        dispatcher.dispatchCommands (commands, 0);

        /** We read the file. */
        u_int64_t checksum = 0;
        size_t nbItemsRead = 0;
        IteratorFile<u_int32_t> it (filename);
        for (it.first(); !it.isDone(); it.next(), nbItemsRead++)  {  checksum += *it;  }

        /** We check the sum of items. */
        CPPUNIT_ASSERT (nbCores*(nbItersPerThread*(nbItersPerThread+1))/2 == checksum);

        //cout << "nbCores=" << nbCores << "  nbItersTotal=" << nbItersTotal << "  cacheSize=" << cacheSize << endl;

        /** We check we read the correct number of items. */
        CPPUNIT_ASSERT (nbItemsRead == nbItersTotal);

        /** We remove the temporary file. */
        System::file().remove (filename);
    }

    /** */
    void bag_checkCache ()
    {
        /** WARNING: nbItersTotal has to be a multiple of nbCores. */

        /** We set the total number of iterations to be done. */
        size_t nbItersTotalTable[] = { 5*1000, 500*1000, 5*1000*1000};

        size_t nbCoresTable[] = { 1, 2, 4, 8};

        /** We set the size of the cache. */
        size_t cacheSizeTable[] = { 10, 1000, 100*1000};

        for (size_t h=0; h<ARRAY_SIZE(nbItersTotalTable); h++)
        {
            for (size_t i=0; i<ARRAY_SIZE(nbCoresTable); i++)
            {
                for (size_t j=0; j<ARRAY_SIZE(cacheSizeTable); j++)
                {
                    bag_checkCache_aux (nbItersTotalTable[h], nbCoresTable[i], cacheSizeTable[j]);
                }
            }
        }
    }
};


/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestBag);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestBag);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
