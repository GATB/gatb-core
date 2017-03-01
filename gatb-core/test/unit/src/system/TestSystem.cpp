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

#include <gatb/system/api/Exception.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/system/impl/MemoryCommon.hpp>
#include <gatb/system/impl/TimeCommon.hpp>
#include <gatb/system/impl/FileSystemCommon.hpp>

#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

#define ABS(a)  ((a)<0 ? -(a) : (a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/********************************************************************************/

/** \brief Test class for operating system operations
 */
class TestSystem : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestSystem);

        CPPUNIT_TEST_GATB (types_check);

        CPPUNIT_TEST_GATB (info_get);

#ifdef WITH_MEMORY_SIZE_STORE
        CPPUNIT_TEST_GATB (memory_basicAlloc);
        CPPUNIT_TEST_GATB (memory_getMemUsage);
        CPPUNIT_TEST_GATB (memory_hugeAlloc);
        CPPUNIT_TEST_GATB (memory_realloc);
        CPPUNIT_TEST_GATB (memory_boundedAllocator);
        CPPUNIT_TEST_GATB (memory_perfAllocator);
#endif
        CPPUNIT_TEST_GATB (memory_memset);
        CPPUNIT_TEST_GATB (memory_memcpy);
        CPPUNIT_TEST_GATB (memory_memcmp);
        // CPPUNIT_TEST_GATB (memory_allocateAll);

        CPPUNIT_TEST_GATB (time_checkSensibility);
        CPPUNIT_TEST_GATB (time_checkException);
        // CPPUNIT_TEST_GATB (time_clockFrequency);  TSC values may be wrong for some arch

        CPPUNIT_TEST_GATB (thread_checkTime);
        CPPUNIT_TEST_GATB (thread_checkSynchro);
        CPPUNIT_TEST_GATB (thread_exception);

        CPPUNIT_TEST_GATB (filesystem_info);
        CPPUNIT_TEST_GATB (filesystem_create_delete);
        CPPUNIT_TEST_GATB (filesystem_iterate);
        CPPUNIT_TEST_GATB (file_create);

//        CPPUNIT_TEST_GATB (file_attributes);

    CPPUNIT_TEST_SUITE_GATB_END();

public:

    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/
    /** \brief check the size of typedef types (u_int32_t for instance)
     *
     * Check types size with sizeof
     */
    void types_check ()
    {
        CPPUNIT_ASSERT (sizeof(u_int8_t)  == 1);
        CPPUNIT_ASSERT (sizeof(u_int16_t) == 2);
        CPPUNIT_ASSERT (sizeof(u_int32_t) == 4);
        CPPUNIT_ASSERT (sizeof(u_int64_t) == 8);

        CPPUNIT_ASSERT (sizeof(int8_t)  == 1);
        CPPUNIT_ASSERT (sizeof(int16_t) == 2);
        CPPUNIT_ASSERT (sizeof(int32_t) == 4);
        CPPUNIT_ASSERT (sizeof(int64_t) == 8);
    }

    /********************************************************************************/
    /** \brief Check that we can retrieve system information.
     *
     * Test of \ref gatb::core::system::ISystemInfo::getNbCores()     \n
     * Test of \ref gatb::core::system::ISystemInfo::getHostName()    \n
     */
    void info_get ()
    {
        CPPUNIT_ASSERT (System::info().getNbCores() >= 1);
        CPPUNIT_ASSERT (System::info().getHostName().empty() == false);
        CPPUNIT_ASSERT (System::info().getVersion().empty() == false);
    }

    /********************************************************************************/
    class Check
    {
    public:

        static Check& singleton()  { static Check instance; return instance; }

        /********************************************************************************/
        void memory_basicAlloc (IMemory& mem)
        {
            /** We snapshot the current memory usage and block numbers. */
            u_int32_t m0 = mem.getCurrentUsage();
            u_int32_t n0 = mem.getNbBlocks();

            void* p1=0;
            void* p2=0;

            /** We check that we have not increased memory usage. */
            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 == 0);
            
            CPPUNIT_ASSERT_NO_THROW (p1 = mem.malloc (    1000));
            CPPUNIT_ASSERT_NO_THROW (p2 = mem.calloc (10, 1000));

            /** We check that we have increased memory usage. */
            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 > 0);

            /** We check that we have two allocated blocks. */
            CPPUNIT_ASSERT (mem.getNbBlocks() - n0 == 2);

            mem.free (p1);
            mem.free (p2);

            /** We check that the memory usage is back as it initial value. */
            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 == 0);

            /** We check that we have no more allocated blocks. */
            CPPUNIT_ASSERT (mem.getNbBlocks() - n0  == 0);
        }

        /********************************************************************************/
        void memory_getMemUsage_aux (IMemory& mem, size_t nbBlocks, size_t* blockSizeTable, size_t sizeTable)
        {
            for (size_t i=0; i<sizeTable; i++)
            {
                size_t blockSize = blockSizeTable[i];

                /** We need a table of pointers. */
                void* ptrs[nbBlocks];

                /** We snapshot the current memory usage. */
                u_int32_t m0 = mem.getCurrentUsage();

                /** We allocate nbBlocks memory blocks. */
                for (size_t i=0; i<nbBlocks; i++)
                {
                    /** We allocate (and keep a reference) a memory block. */
                    ptrs [i] =  mem.malloc (blockSize);

                    //printf ("+++ [%d]  %d  %d\n", i, mem.getCurrentUsage() - m0,  (i+1)*blockSize);

                    /** We check that the current memory usage is increasing. */
                    CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 >= (i+1)*blockSize);
                }

                for (size_t i=0; i<nbBlocks; i++)
                {
                    /** We deallocate the memory block. */
                    mem.free (ptrs[i]);

                    //printf ("--- [%d]  %d  %d\n", i, mem.getCurrentUsage() - m0,  (nbBlocks-i-1)*blockSize);

                    /** We check that the current memory usage is decreasing. */
                    CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 < (nbBlocks-i)*blockSize );
                }

                /** We check that the current memory usage is the same at the beginning of the test. */
                CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 == 0);
            }
        }

        /********************************************************************************/
        void memory_hugeAlloc (IMemory& mem)
        {
            void* p1=0;
            void* p2=0;

            /** We snapshot the current memory usage and block numbers. */
            u_int32_t m0 = mem.getCurrentUsage();
            u_int32_t n0 = mem.getNbBlocks();

            CPPUNIT_ASSERT_NO_THROW (p1 = mem.malloc (KBYTE+0));
            CPPUNIT_ASSERT_THROW    (p2 = mem.malloc (KBYTE+1), gatb::core::system::Exception);

            CPPUNIT_ASSERT (mem.getNbBlocks() - n0  == 1);

            CPPUNIT_ASSERT (p1 != 0);
            CPPUNIT_ASSERT (p2 == 0);

            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 > 0);
            mem.free (p1);
            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 == 0);

            CPPUNIT_ASSERT (mem.getNbBlocks() - n0  == 0);
        }

        /********************************************************************************/
        void memory_realloc (IMemory& mem)
        {
            void* p1=0;

            CPPUNIT_ASSERT (p1 == 0);

            for (size_t i=1; i<=100; i++)
            {
                p1 = mem.realloc (p1, i*1000);
                CPPUNIT_ASSERT (p1 != 0);
            }

            mem.free (p1);
        }

        /********************************************************************************/
        void memory_boundedAllocator (IMemory& mem)
        {
            /** We snapshot the current memory usage and block numbers. */
            u_int32_t m0 = mem.getCurrentUsage();
            u_int32_t n0 = mem.getNbBlocks();

            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 == 0);
            CPPUNIT_ASSERT (mem.getNbBlocks()     - n0 == 0);

            void* p1 = mem.malloc (130);
            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 >= 130);
            CPPUNIT_ASSERT (mem.getNbBlocks()     - n0 == 1);

            mem.free (p1);
            CPPUNIT_ASSERT (mem.getCurrentUsage() - m0 == 0);
            CPPUNIT_ASSERT (mem.getNbBlocks()     - n0 == 0);
        }

        /********************************************************************************/
        ITime::Value memory_perfAllocator (IMemoryAllocator& mem, size_t nb, IMemory::BlockSize_t size)
        {
            TimeSystem usec (ITime::MSEC);

            ITime::Value t0 = usec.getTimeStamp();

            for (size_t i=0; i<nb; i++)
            {
                void* ptr = mem.malloc (size);
                mem.free (ptr);
            }

            ITime::Value t1 = usec.getTimeStamp();

            return t1 - t0;
        }
    };

    /********************************************************************************/
    /********************************************************************************/
    /********************************************************************************/
    /** \brief Check basic allocations and check statistics accordingly
     *
     *  Test of \ref gatb::core::system::IMemory::getCurrentUsage()   \n
     *  Test of \ref gatb::core::system::IMemory::getNbBlocks()       \n
     *  Test of \ref gatb::core::system::IMemory::malloc()            \n
     *  Test of \ref gatb::core::system::IMemory::calloc()            \n
     *  Test of \ref gatb::core::system::IMemory::free()              \n
     */
    void memory_basicAlloc ()
    {

        MemoryBounded allocator (System::memory(), 2048*KBYTE, 11*2048*KBYTE);
        MemoryCommon mem (allocator, MemoryOperationsCommon::singleton());

        Check::singleton().memory_basicAlloc (System::memory());
        Check::singleton().memory_basicAlloc (mem);
    }

    /********************************************************************************/
    /** \brief Test of memory usage evolution
     *
     *  1) allocate several blocks and check that the used memory increases     \n
     *  2) deallocate the blocks and check that the used memory decreases       \n
     *  3) check that the final memory usage is the same as at the test start   \n
     *
     * Test of \ref gatb::core::system::IMemory::getCurrentUsage()    \n
     * Test of \ref gatb::core::system::IMemory::malloc()             \n
     * Test of \ref gatb::core::system::IMemory::free()               \n
     */
    void memory_getMemUsage ()
    {
        size_t nbBlocks = 10;
        size_t blockSizeTable[] = { 128*KBYTE, 1024*KBYTE, 2048*KBYTE };
        size_t blockSizeLength  = sizeof(blockSizeTable) / sizeof(blockSizeTable[0]);

        /** We create a bounded memory allocator. */
        MemoryBounded allocator (System::memory(), 2048*KBYTE, 11*2048*KBYTE);
        MemoryCommon mem (allocator, MemoryOperationsCommon::singleton());

        /** Note: since the getMemUsage is not completely accurate for the generic allocator,
         * we use big block sizes in order not to have failed tests here. */

        Check::singleton().memory_getMemUsage_aux (System::memory(), nbBlocks, blockSizeTable, blockSizeLength);
        Check::singleton().memory_getMemUsage_aux (mem,                nbBlocks, blockSizeTable, blockSizeLength);
    }

    /********************************************************************************/
    /** \brief Force to have an exception due to too big required size
     *
     * This test checks that an impossible malloc raises an exception.
     *
     * Test of \ref gatb::core::system::IMemory::getCurrentUsage()    \n
     * Test of \ref gatb::core::system::IMemory::getNbBlocks()        \n
     * Test of \ref gatb::core::system::IMemory::malloc()             \n
     * Test of \ref gatb::core::system::IMemory::free()               \n
     */
    void memory_hugeAlloc ()
    {
        MemoryBounded allocator (System::memory(), KBYTE, 100*KBYTE);
        MemoryCommon mem (allocator, MemoryOperationsCommon::singleton());

        Check::singleton().memory_hugeAlloc (mem);
    }

    /********************************************************************************/
    /** \brief Reallocation test
     *
     * Test of \ref gatb::core::system::IMemory::realloc()    \n
     * Test of \ref gatb::core::system::IMemory::free()       \n
     */
    void memory_realloc ()
    {
        MemoryBounded allocator (System::memory(), MBYTE, 100*MBYTE);
        MemoryCommon mem (allocator, MemoryOperationsCommon::singleton());

        Check::singleton().memory_realloc (System::memory());
        Check::singleton().memory_realloc (mem);
    }

    /********************************************************************************/
    /** \brief Test of the \ref gatb::core::system::impl::MemoryBounded class
     *
     *  Test of \ref gatb::core::system::IMemory::getCurrentUsage()   \n
     *  Test of \ref gatb::core::system::IMemory::getNbBlocks()       \n
     *  Test of \ref gatb::core::system::IMemory::malloc()            \n
     *  Test of \ref gatb::core::system::IMemory::free()              \n
     */
    void memory_boundedAllocator ()
    {
        MemoryBounded allocator (System::memory(), 128*KBYTE, 64*MBYTE);
        MemoryCommon mem (allocator, MemoryOperationsCommon::singleton());

        Check::singleton().memory_boundedAllocator (mem);
    }

    /********************************************************************************/
    /** \brief Memory allocation performance test
     *
     * Test of \ref gatb::core::system::IMemory::malloc() \n
     * Test of \ref gatb::core::system::IMemory::free()   \n
     * Test of \ref gatb::core::system::impl::TimeSystem::getTimeStamp()   \n
     */
    void memory_perfAllocator ()
    {
        size_t               nb   = 10*1000;
        IMemory::BlockSize_t size = 64*MBYTE;

        Check::singleton().memory_perfAllocator (System::memory(),                   nb, size);
        Check::singleton().memory_perfAllocator (MemoryAllocatorStdlib::singleton(), nb, size);
    }

    /********************************************************************************/
    /** \brief Check memset operation
     *
     * Test of \ref gatb::core::system::IMemory::malloc() \n
     * Test of \ref gatb::core::system::IMemory::memset() \n
     * Test of \ref gatb::core::system::IMemory::free()   \n
     */
    void memory_memset ()
    {
        size_t nb = 10000;
        int c = 47;

        /** We allocate a block. */
        u_int8_t* ptr = (u_int8_t*) System::memory().malloc (nb);
        CPPUNIT_ASSERT (ptr != 0);

        /** We reset all values of the block. */
        System::memory().memset (ptr, c, nb);

        /** We check that all values are ok. */
        for (size_t i=0; i<nb; i++)  {  CPPUNIT_ASSERT (ptr[i] == c); }

        /** We release the block. */
        System::memory().free (ptr);
    }

    /********************************************************************************/
    /** \brief Check memcpy operation
     *
     * Test of \ref gatb::core::system::IMemory::malloc() \n
     * Test of \ref gatb::core::system::IMemory::memset() \n
     * Test of \ref gatb::core::system::IMemory::memcpy() \n
     * Test of \ref gatb::core::system::IMemory::free()   \n
     */
    void memory_memcpy ()
    {
        size_t nb = 10000;
        int c1=11, c2=47;

        /** We allocate a block. */
        u_int8_t* ptr1 = (u_int8_t*) System::memory().malloc (nb);
        CPPUNIT_ASSERT (ptr1 != 0);

        /** We allocate a block. */
        u_int8_t* ptr2 = (u_int8_t*) System::memory().malloc (nb);
        CPPUNIT_ASSERT (ptr2 != 0);

        /** We reset all values of the blocks. */
        System::memory().memset (ptr1, c1, nb);
        System::memory().memset (ptr2, c2, nb);

        /** We check that all values are ok. */
        for (size_t i=0; i<nb; i++)  {  CPPUNIT_ASSERT (ptr1[i] == c1); }
        for (size_t i=0; i<nb; i++)  {  CPPUNIT_ASSERT (ptr2[i] == c2); }

        /** We copy block 1 into block 2. */
        System::memory().memcpy (ptr2, ptr1, nb);

        /** We check that all values are ok. */
        for (size_t i=0; i<nb; i++)  {  CPPUNIT_ASSERT (ptr2[i] == c1); }

        /** We release the blocks. */
        System::memory().free (ptr1);
        System::memory().free (ptr2);
    }

    /********************************************************************************/
    /** \brief Check memcmp operation
     *
     * Test of \ref gatb::core::system::IMemory::malloc() \n
     * Test of \ref gatb::core::system::IMemory::memcmp() \n
     * Test of \ref gatb::core::system::IMemory::free()   \n
     */
    void memory_memcmp ()
    {
        size_t nb = 10000;
        int c1=11, c2=47;

        /** We allocate a block. */
        u_int8_t* ptr1 = (u_int8_t*) System::memory().malloc (nb);
        CPPUNIT_ASSERT (ptr1 != 0);

        /** We allocate a block. */
        u_int8_t* ptr2 = (u_int8_t*) System::memory().malloc (nb);
        CPPUNIT_ASSERT (ptr2 != 0);

        /** We reset all values of the blocks. */
        System::memory().memset (ptr1, c1, nb);
        System::memory().memset (ptr2, c2, nb);

        CPPUNIT_ASSERT (System::memory().memcmp (ptr1, ptr2, nb) != 0);

        /** We copy block 1 into block 2. */
        System::memory().memcpy (ptr2, ptr1, nb);

        /** We check that all values are ok. */
        CPPUNIT_ASSERT (System::memory().memcmp (ptr1, ptr2, nb) == 0);

        /** We release the blocks. */
        System::memory().free (ptr1);
        System::memory().free (ptr2);
    }

    /********************************************************************************/
    void memory_allocateAll ()
    {
        list<void*> ptrs;

        IMemory::BlockSize_t len   =  1ull*MBYTE;
        IMemory::TotalSize_t total =  6ull*GBYTE;

        size_t nb = 0;

        MemoryBounded allocator (System::memory(), 1*len, total);
        MemoryCommon mem (allocator, MemoryOperationsCommon::singleton());

        while (true)
        {
            try
            {
                void* ptr = mem.malloc (len);

                if (ptr != 0)
                {
                    mem.memset (ptr, 0, len);

                    ptrs.push_back (ptr);
                    nb++;
                }
            }
            catch (...)
            {
                break;
            }
        }

        for (list<void*>::iterator it = ptrs.begin(); it != ptrs.end(); it++)
        {
            mem.free (*it);
        }
    }

    /********************************************************************************/
    /** \brief Check that ITime sensibility
     *
     *  This test waits for some time and checks that the elapsed time is coherent with the
     *  retrieved values from several TimeSystem instances.
     *
     *  Test of \ref gatb::core::system::ITime::getTimeStamp()   \n
     */
    void time_checkSensibility ()
    {
        TimeSystem usec (ITime::USEC);
        TimeSystem msec (ITime::MSEC);
        TimeSystem  sec (ITime::SEC);

        list<ITime*> timeList;
        timeList.push_back (&usec);
        timeList.push_back (&msec);
        timeList.push_back (&sec);

        size_t delayInSecond = 5;

        for (list<ITime*>::iterator it = timeList.begin(); it != timeList.end(); it++)
        {
            ITime::Value t0 = (*it)->getTimeStamp();
            sleep (delayInSecond);
            ITime::Value t1 = (*it)->getTimeStamp();

            CPPUNIT_ASSERT ((t1-t0) / (*it)->getUnit() == delayInSecond);
        }
    }

    /********************************************************************************/
    /** \brief Check that class TimeSystem is correcty constructed.
     *
     * Test of \ref gatb::core::system::impl::TimeSystem  \n
     */
    void time_checkException ()
    {
        CPPUNIT_ASSERT_NO_THROW (TimeSystem s (ITime::USEC));
        CPPUNIT_ASSERT_NO_THROW (TimeSystem s (ITime::MSEC));
        CPPUNIT_ASSERT_NO_THROW (TimeSystem s (ITime::SEC));
        CPPUNIT_ASSERT_THROW    (TimeSystem s (ITime::UNDEFINED), gatb::core::system::Exception);
    }

    /********************************************************************************/
    /** \brief Check that class TimeCycle retrieves correctly
     * an estimation of the CPU frequency clock.
     *
     * Test of \ref gatb::core::system::impl::TimeCycle::getClockFrequency()  \n
     */
    void time_clockFrequency ()
    {
        /** We create a ITime instance with CPU cycle measures. */
        TimeCycle ts;

        CPPUNIT_ASSERT (ts.getClockFrequency() > 0);
    }

    /********************************************************************************/
    static void* thread_checkTime_mainloop (void* data)
    {
        ITime::Value* totalTime = (ITime::Value*)data;

        ITime::Value t0 = System::time().getTimeStamp();
        sleep (2);
        ITime::Value t1 = System::time().getTimeStamp();

        __sync_fetch_and_add (totalTime, t1 - t0);

        return 0;
    }

    /********************************************************************************/
    /** \brief Check threads behavior
     *
     *  This test creates N threads and each thread waits for T seconds.
     *  At the end, the aggregated time spent in the threads should be N.T and the
     *  user time should be T. Note that the check is done by computing
     *  (threadsTime - nbThreads*userTime) with some tolerance (up do 5%) because
     *  it is unlikely that the relative error will always be 0.
     *
     *  The test is done each time with 1, 2, ... up to the number of available cores.
     *
     *  Test of \ref gatb::core::system::ISystemInfo::getNbCores()    \n
     *  Test of \ref gatb::core::system::IThreadFactory::newThread()  \n
     *  Test of \ref gatb::core::system::IThread::join()              \n
     *  Test of \ref gatb::core::system::ITime::getTimeStamp()        \n
     */
    void thread_checkTime ()
    {
        /** We loop over the interval [1..nbCores] */
        for (size_t nbCores=1; nbCores<=System::info().getNbCores(); nbCores++)
        {
            ITime::Value totalTime = 0;

            /** We create some threads. */
            IThread* threads [nbCores];
            for (size_t i=0; i<nbCores; i++)   {  threads[i] = System::thread().newThread (thread_checkTime_mainloop, &totalTime);  }

            ITime::Value t0 = System::time().getTimeStamp();

            /** We wait the end of each thread. */
            for (size_t i=0; i<nbCores; i++)   {  threads[i]->join();    delete threads[i];  }

            ITime::Value t1 = System::time().getTimeStamp();

            /** We compute the difference between threads time and user time. */
            double err = 100.0 * ABS (1.0 - (double) (nbCores*(t1-t0)) / (double) totalTime);

            //printf ("[%d]  %d  %d  %d  err=%lf \n", nbCores, totalTime, t1-t0, totalTime - nbCores*(t1-t0), err);

            /** We check that the relative error is not too big. Note that 5% here allows to run the test even if
             * the CPU is already widely used by other processes; since the threads main loop only just sleeps for
             * a while (and is not time consuming), this 5% margin should be enough not to have false negative test
             * report. */
            CPPUNIT_ASSERT (err < 5.0);
        }
    }

    /********************************************************************************/

    struct Data
    {
        Data (size_t nb, ISynchronizer* s=0) : nbIter(nb), value(0), synchro(s)  {}

        size_t         nbIter;
        size_t         value;
        ISynchronizer* synchro;
    };

    static void* thread_checkSynchro_mainloop (void* arg)
    {
        Data& data = *(Data*) arg;

        for (size_t i=0; i<data.nbIter; i++)
        {
            if (data.synchro)  {  data.synchro->lock(); }

            /** We compute the new value. */
            data.value ++;

            /** We simulate the fact that each thread takes different time for the computation. */
            usleep (rand() % 10 + 1);

            if (data.synchro)  {  data.synchro->unlock(); }
        }

        return 0;
    }

    /********************************************************************************/
    /** \brief Check threads and synchronizers behavior
     *
     *  This test creates N threads and each thread increase C times a counter.
     *  At the end, the counter value should be N.C
     *
     *  In order to make the test work, we need to use lock/unlock in each thread, otherwise
     *  the final counter value would be wrong.
     *
     *  Test of \ref gatb::core::system::ISystemInfo::getNbCores()            \n
     *  Test of \ref gatb::core::system::IThreadFactory::newThread()          \n
     *  Test of \ref gatb::core::system::IThreadFactory::newSynchronizer()    \n
     *  Test of \ref gatb::core::system::IThread::join()                      \n
     */
    void thread_checkSynchro ()
    {
        /** We create a synchronizer to be shared by the different threads. */
        ISynchronizer* synchro = System::thread().newSynchronizer();

        /** We initialize our global resource. */
        Data data (1000, synchro);

        /** We retrieve the number of cores. */
        size_t nbCores = System::info().getNbCores();
        CPPUNIT_ASSERT (nbCores > 0);

        /** We create some threads. */
        IThread* threads [nbCores];
        for (size_t i=0; i<nbCores; i++)   {  threads[i] = System::thread().newThread (thread_checkSynchro_mainloop, &data);  }

        /** We wait the end of each thread. */
        for (size_t i=0; i<nbCores; i++)   {  threads[i]->join();    delete threads[i];  }

        CPPUNIT_ASSERT (data.value == nbCores * data.nbIter);

        /** Some cleanup. */
        delete synchro;
    }

    /********************************************************************************/
    static void* thread_exception_mainloop (void* arg)
    {
        IThreadGroup* threadGroup = (IThreadGroup*) arg;

        /** We have to catch the thread local exception and forward it to the thread group. */
        try
        {
            /** We launch an exception. */
            throw core::system::Exception ("something wrong");
        }
        catch (core::system::Exception& e)
        {
            threadGroup->addException (e);
        }

        return 0;
    }

    void thread_exception ()
    {
        /** We create a thread group. */
        IThreadGroup* threadGroup = ThreadGroup::create ();

        /** We retrieve the number of cores. */
        size_t nbCores = System::info().getNbCores();
        CPPUNIT_ASSERT (nbCores > 0);

        /** We add some threads to the group. */
        for (size_t i=0; i<nbCores; i++)   {  threadGroup->add (thread_exception_mainloop, threadGroup);  }

        /** We start the group. */
        threadGroup->start ();

        /** Now, the treads are all finished and joined, we check for got exceptions. */
        CPPUNIT_ASSERT (threadGroup->hasExceptions());

        try
        {
            throw threadGroup->getException();
            CPPUNIT_ASSERT (false);
        }
        catch (core::system::Exception& e)
        {
            CPPUNIT_ASSERT (true);
        }

        /** We have to destroy the thread group this way. */
        ThreadGroup::destroy(threadGroup);
    }

    /********************************************************************************/
    /** \brief check information from the file system.
     *
     *  Test of \ref gatb::core::system::IFileSystem::getMaxFilesNumber()     \n
     *  Test of \ref gatb::core::system::IFileSystem::getCurrentDirectory()   \n
     *  Test of \ref gatb::core::system::IFileSystem::getTemporaryDirectory() \n
     *  Test of \ref gatb::core::system::IFileSystem::doesExist()             \n
     *  Test of \ref gatb::core::system::IFileSystem::getAvailableSpace()     \n
     */
    void filesystem_info ()
    {
        CPPUNIT_ASSERT (System::file().getMaxFilesNumber() > 0);

        CPPUNIT_ASSERT (System::file().getCurrentDirectory().empty() == false);
        CPPUNIT_ASSERT (System::file().doesExist (System::file().getCurrentDirectory()));

        CPPUNIT_ASSERT (System::file().getTemporaryDirectory().empty() == false);
        CPPUNIT_ASSERT (System::file().doesExist (System::file().getTemporaryDirectory()));

        CPPUNIT_ASSERT (System::file().getAvailableSpace (System::file().getCurrentDirectory()) > 0);
    }

    /********************************************************************************/
    /** \brief creation and deletion of directory
     *
     * We create a directory in the temporary directory; its name is the current date.
     *
     * Test of \ref gatb::core::system::IFileSystem::getTemporaryDirectory()  \n
     * Test of \ref gatb::core::system::IFileSystem::mkdir()                  \n
     * Test of \ref gatb::core::system::IFileSystem::doesExist()              \n
     * Test of \ref gatb::core::system::IFileSystem::rmdir()                  \n
     * Test of \ref gatb::core::system::ITime::getDateString()                \n
     *
     */
    void filesystem_create_delete ()
    {
        /** We retrieve the temporary directory. */
        string tmpdir = System::file().getTemporaryDirectory();
        CPPUNIT_ASSERT (tmpdir.empty() == false);

        /** We build a test directory path. */
        stringstream ss;   ss << tmpdir << "/" << System::time().getDateString();
        IFileSystem::Path dirpath = ss.str();

        /** We create a directory. */
        int res1 = System::file().mkdir (dirpath, 01777);
        CPPUNIT_ASSERT (res1 == 0);

        /** We check that the directory exists. */
        CPPUNIT_ASSERT (System::file().doesExist (dirpath) == true);

        /** We delete the directory. */
        int res2 = System::file().rmdir (dirpath);
        CPPUNIT_ASSERT (res2 == 0);

        /** We check that the directory exists. */
        CPPUNIT_ASSERT (System::file().doesExist (dirpath) == false);
    }

    /********************************************************************************/
    static void iterateAction (const IFileSystem::Path& entry, void* data)
    {
        size_t& nb = *(size_t*)data;
        nb++;
    }

    /** \brief iterate some directory content
     *
     * Test of \ref gatb::core::system::IFileSystem::getCurrentDirectory()    \n
     * Test of \ref gatb::core::system::IFileSystem::iterate()                \n
     */
    void filesystem_iterate ()
    {
        size_t nbEntries = 0;

        /** We retrieve the current directory. */
        string currentdir = System::file().getCurrentDirectory ();
        CPPUNIT_ASSERT (currentdir.empty() == false);

        /** We iterate the directory. */
        System::file().iterate (currentdir, iterateAction, &nbEntries);

        CPPUNIT_ASSERT (nbEntries > 0);
    }

    /********************************************************************************/
    /** \brief Creation of a file with some content
     *
     * Test of \ref gatb::core::system::IFileSystem::getCurrentDirectory() \n
     * Test of \ref gatb::core::system::IFileSystem::newFile()             \n
     * Test of \ref gatb::core::system::IFile::isOpen()   \n
     * Test of \ref gatb::core::system::IFile::print()    \n
     * Test of \ref gatb::core::system::IFile::seeko()    \n
     * Test of \ref gatb::core::system::IFile::tell()     \n
     * Test of \ref gatb::core::system::IFile::get()      \n
     * Test of \ref gatb::core::system::IFile::isEOF()    \n
     */
    void file_create ()
    {
        const char* info = "0123456789abcdefghijklmnopqrstuvwxyz";

        /** We retrieve the temporary directory. */
        string tmpdir = System::file().getTemporaryDirectory();
        CPPUNIT_ASSERT (tmpdir.empty() == false);

        /** We create the file handle. We need it first for write. */
        IFile* writer = System::file().newFile (tmpdir, "dummy", "w");

        /** We check the file handle exists. Note that we will ALWAYS get an instance,
         * even if the path is incorrect for instance. */
        CPPUNIT_ASSERT (writer != 0);

        if (writer->isOpen())
        {
            /** We add some information into the file. */
            writer->print (info);

            /** We check the size of the file. */
            CPPUNIT_ASSERT (writer->getSize() == strlen (info));
        }

        /** We create the file handle. We need it for read. */
        IFile* reader = System::file().newFile (writer->getPath(), "r");

        /** We check the file handle exists. Note that we will ALWAYS get an instance,
         * even if the path is incorrect for instance. */
        CPPUNIT_ASSERT (reader != 0);

        if (reader->isOpen())
        {
            /** We are supposed to be at the beginning of the file. */
            CPPUNIT_ASSERT (reader->tell() == 0);

            /** We get one character and check it is ok. */
            const char* loop = info;
            for (char c=reader->get();  ! reader->isEOF();  c=reader->get(), loop++)
            {
                CPPUNIT_ASSERT (c == *loop);
            }

            CPPUNIT_ASSERT (reader->tell() == strlen (info));
        }

        /** We check that the file exists */
        CPPUNIT_ASSERT (System::file().doesExist (writer->getPath()));

        /** We remove the file. */
        System::file().remove (writer->getPath());

        /** We check that the file doesn't exist any more */
        CPPUNIT_ASSERT (System::file().doesExist (writer->getPath()) == false);

        /** Some cleanup. */
        delete writer;
        delete reader;
    }

    /********************************************************************************/
    void file_attributes ()
    {
        /** We create a temporary file. */
        string filename = System::file().getTemporaryDirectory() + "/foo.txt";

        IFile* file = System::file().newFile (filename, "w");
        CPPUNIT_ASSERT (file != 0);
        CPPUNIT_ASSERT (System::file().doesExist (filename) == true);

        const char* key = "myKey";
        int res = 0;

        /** We set some tags to the file. */
        res = System::file().setAttribute (filename, key, "%d", 147);
        CPPUNIT_ASSERT (res >= 0);

        /** We get the attribute. */
        string value;
        res = System::file().getAttribute (filename, key, value);

        CPPUNIT_ASSERT (res == 3);

        /** We delete the temporary file. */
        System::file().remove (filename);

        CPPUNIT_ASSERT (System::file().doesExist (filename) == false);

        /** Cleanup */
        delete file;
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestSystem);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestSystem);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/

