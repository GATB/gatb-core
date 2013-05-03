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

#include <gatb/tools/math/ttmath/ttmath.h>
#include <gatb/tools/math/LargeInt.hpp>

using namespace std;
using namespace gatb::core::tools::math;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for miscellaneous operations
 */
class TestMath : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestMath);

        CPPUNIT_TEST_GATB (math_checkBasic);
        CPPUNIT_TEST_GATB (math_checkFibo);

    CPPUNIT_TEST_SUITE_GATB_END();

public:
    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    template <typename T> void math_checkBasicTemplate ()
    {
        T a (3);
        T b (5);

        CPPUNIT_ASSERT (a == 3);
        CPPUNIT_ASSERT (b == 5);
        CPPUNIT_ASSERT (a != b);
        CPPUNIT_ASSERT (a <  b);
        CPPUNIT_ASSERT (a <= b);

        CPPUNIT_ASSERT (a + 1 == 4);
        CPPUNIT_ASSERT (b + 1 == 6);

        CPPUNIT_ASSERT (a - 1 == 2);
        CPPUNIT_ASSERT (b - 1 == 4);

        CPPUNIT_ASSERT (a*2 == 6);
        CPPUNIT_ASSERT (b*4 == 20);

        CPPUNIT_ASSERT (a/2 == 1);
        CPPUNIT_ASSERT (b/2 == 2);

        CPPUNIT_ASSERT (a%2 == 1);
        CPPUNIT_ASSERT (b%2 == 1);

        CPPUNIT_ASSERT ( (a ^ b) == 6);
        CPPUNIT_ASSERT ( (a & b) == 1);

        CPPUNIT_ASSERT ( (a << 4) == 48);
        CPPUNIT_ASSERT ( (a >> 1) ==  1);
        CPPUNIT_ASSERT ( (a >> 2) ==  0);
    }

    /********************************************************************************/
    template <typename T> void math_checkFiboTemplate ()
    {
        static long table[] =
        {
            1, 2, 3, 5, 8, 13, 21, 34, 55,
            89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765,
            10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040,
            1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155,
            165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073, 4807526976, 7778742049, 12586269025
        };

        T u0(0), u1(1);

        for (size_t i=0; i<sizeof(table)/sizeof(table[0]); i++)
        {
            T u2;
            u2 = u1 + u0;
            CPPUNIT_ASSERT (u2 == table[i]);

            u0 = u1;
            u1 = u2;
        }
    }

    /********************************************************************************/
    void math_checkBasic ()
    {
        math_checkBasicTemplate < u_int64_t >                    ();
        math_checkBasicTemplate < ttmath::UInt<KMER_PRECISION> > ();
        math_checkBasicTemplate < LargeInt<KMER_PRECISION> >     ();
    }

    /********************************************************************************/
    void math_checkFibo ()
    {
        math_checkFiboTemplate < u_int64_t >                    ();
        math_checkFiboTemplate < ttmath::UInt<KMER_PRECISION> > ();
        math_checkFiboTemplate < LargeInt<KMER_PRECISION> >     ();
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION (TestMath);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
