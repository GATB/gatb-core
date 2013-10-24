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

#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/tools/math/Integer.hpp>

using namespace std;
using namespace gatb::core::tools::math;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

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
    template <typename T, typename U=T> void math_checkBasicTemplate ()
    {
#define CST(a) T(U(a))

        T a (CST(3));
        T b (CST(5));

        CPPUNIT_ASSERT (a == CST(3));
        CPPUNIT_ASSERT (b == CST(5));

        CPPUNIT_ASSERT (a != b);
        CPPUNIT_ASSERT (a <  b);
        CPPUNIT_ASSERT (a <= b);

        CPPUNIT_ASSERT (a + CST(1) == CST(4));
        CPPUNIT_ASSERT (b + CST(1) == CST(6));

        CPPUNIT_ASSERT (a - CST(1) == CST(2));
        CPPUNIT_ASSERT (b - CST(1) == CST(4));

        CPPUNIT_ASSERT (a*2 == CST(6));
        CPPUNIT_ASSERT (b*4 == CST(20));

        CPPUNIT_ASSERT (a/2 == CST(1));
        CPPUNIT_ASSERT (b/2 == CST(2));

        CPPUNIT_ASSERT (a%2 == 1);
        CPPUNIT_ASSERT (b%2 == 1);

        CPPUNIT_ASSERT ( (a ^ b) == CST(6));
        CPPUNIT_ASSERT ( (a & b) == CST(1));

        CPPUNIT_ASSERT ( (a << 4) == CST(48));
        CPPUNIT_ASSERT ( (a >> 1) == CST(1));
        CPPUNIT_ASSERT ( (a >> 2) == CST(0));

        a  = CST(4);
        a += CST(5);
        CPPUNIT_ASSERT (a == CST(9));

        a = CST (0xCC);
        b = ~a  & CST(0xFF);
        CPPUNIT_ASSERT (b == CST(0x33));

        CPPUNIT_ASSERT (revcomp (CST (0x112233445566), 11) ==  CST (0xcffee));
#if 0
        CPPUNIT_ASSERT (simplehash16 (CST(0x11223344), 12) == 2129748964359514508);
        CPPUNIT_ASSERT (hash1(CST(0x11223344), 5)  ==  9479217074923583263UL);   // => BUG ???
        CPPUNIT_ASSERT (oahash(CST(0x11223344))    == 102913755990806762);       // => BUG ???
#endif
    }

    /********************************************************************************/
    template <typename T, typename U=T> void math_checkFiboTemplate ()
    {
        static u_int64_t table[] =
        {
            1, 2, 3, 5, 8, 13, 21, 34, 55,
            89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765,
            10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040,
            1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155,
            165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073, 4807526976u, 7778742049u, 12586269025u
        };

        T u0(0), u1(1);

        for (size_t i=0; i<sizeof(table)/sizeof(table[0]); i++)
        {
            T u2(0);
            u2 = u1 + u0;
            CPPUNIT_ASSERT (u2 == T(table[i]));

            u0 = u1;
            u1 = u2;
        }
    }

    /********************************************************************************/
    void math_checkBasic ()
    {
        math_checkBasicTemplate < LargeInt<1> >     ();
        math_checkBasicTemplate < LargeInt<2> >     ();
        math_checkBasicTemplate < LargeInt<3> >     ();
        math_checkBasicTemplate < LargeInt<4> >     ();
        math_checkBasicTemplate < LargeInt<5> >     ();
        math_checkBasicTemplate < LargeInt<6> >     ();

        math_checkBasicTemplate <Integer, LargeInt<1> > ();
        math_checkBasicTemplate <Integer, LargeInt<2> > ();
        math_checkBasicTemplate <Integer, LargeInt<3> > ();
        math_checkBasicTemplate <Integer, LargeInt<4> > ();
    }

    /********************************************************************************/
    void math_checkFibo ()
    {
        math_checkBasicTemplate < LargeInt<1> >     ();
        math_checkBasicTemplate < LargeInt<2> >     ();
        math_checkBasicTemplate < LargeInt<3> >     ();
        math_checkBasicTemplate < LargeInt<4> >     ();

        math_checkBasicTemplate <Integer, LargeInt<1> > ();
        math_checkBasicTemplate <Integer, LargeInt<2> > ();
        math_checkBasicTemplate <Integer, LargeInt<3> > ();
        math_checkBasicTemplate <Integer, LargeInt<4> > ();
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestMath);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestMath);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
