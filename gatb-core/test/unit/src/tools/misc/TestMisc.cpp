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
#include <gatb/system/impl/TimeCommon.hpp>

#include <gatb/tools/misc/api/Range.hpp>
#include <gatb/tools/misc/api/Vector.hpp>
#include <gatb/tools/misc/api/Macros.hpp>

#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/misc/impl/StringLine.hpp>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <memory>

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

#define ABS(a)  ((a)<0 ? -(a) : (a))

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for miscellaneous operations
 */
class TestMisc : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestMisc);

        CPPUNIT_TEST_GATB (range_checkIterator1);
        CPPUNIT_TEST_GATB (range_checkIterator2);
        //  NEED A GOOD TIMER ACCURACY...  CPPUNIT_TEST_GATB (range_checkPerformance);
        // DEACTIVATED BECAUSE OF MACOS (TO BE INVESTIGATED...)  CPPUNIT_TEST_GATB (vector_check1);
        CPPUNIT_TEST_GATB (vector_check2);
        CPPUNIT_TEST_GATB (vector_check3);
        CPPUNIT_TEST_GATB (parser_check1);
        CPPUNIT_TEST_GATB (parser_check2);

        CPPUNIT_TEST_GATB (stringline_check1);

    CPPUNIT_TEST_SUITE_GATB_END();

public:
    /********************************************************************************/
    void setUp    ()  {  srand (time(NULL));  }
    void tearDown ()  {}

    /********************************************************************************/
    /** \brief test of Range class and its associated iterator.
     *
     * We create an iterator for an integer range and check that each iterated item is ok.
     */
    void range_checkIterator1 ()
    {
        /** We create a range. */
        Range<size_t> range (1, 1000);

        /** We create an iterator from the range. */
        Range<size_t>::Iterator it (range);

        size_t check = range.getBegin();

        /** We iterate each item of the range. */
        for (it.first(); !it.isDone(); it.next())
        {
            /** We check that the current iterated item is ok. */
            CPPUNIT_ASSERT (it.item() == check++);
        }
    }

    /********************************************************************************/
    /** \brief test of Range associated iterator.
     *
     * We create an iterator for an integer range and check that each iterated item is ok.
     */
    void range_checkIterator2 ()
    {
        size_t from=1, to=100, check=from;

        /** We create an iterator from the range. */
        Range<size_t>::Iterator it (from, to);

        /** We iterate each item of the range. */
        for (it.first(); !it.isDone(); it.next())
        {
            /** We check that the current iterated item is ok. */
            CPPUNIT_ASSERT (it.item() == check++);
        }
    }

    /********************************************************************************/
    void range_checkPerformance ()
    {
        typedef u_int64_t INT;

        u_int64_t sum1=0, sum2=0;

        INT begin=1, end=1000*1000*1000;

        /** We create a range. */
        Range<INT> range (begin, end);

        /** We create an iterator from the range (seen as an Iterable). */
        Range<INT>::Iterator it (range);

        /** We need something to measure elapsed time. */
        ITime& ts = System::time();

        /** Note that we use 'volatile' keyword in order to avoid unwanted optimizations here. */
        volatile ITime::Value t0 = ts.getTimeStamp();

        for (it.first(); !it.isDone(); it.next())   {  sum1 += it.item();  }

        volatile ITime::Value t1 = ts.getTimeStamp();

        for (INT i=begin; i<=end; i++)   {  sum2 += i;  }

        volatile ITime::Value t2 = ts.getTimeStamp();

        /** We check we got the same result with both loops. */
        CPPUNIT_ASSERT (sum1 > 0  &&  sum1 == sum2);

        /** We check that performances are of same order. */
        double err = 100.0 * ABS (1.0 - (double)(t1-t0) / (double)(t2-t1));

        /** The 'iterator' loop should be less than 100% slower than the 'direct' loop. */
        CPPUNIT_ASSERT (err < 100);
    }

    /********************************************************************************/
    void vector_check1 ()
    {
        char table[] = { 1, 2, 3, 5, 8, 13, 21, 34, 55, 89};

        Vector<char> v1;
        CPPUNIT_ASSERT (v1.size() == 0);
        CPPUNIT_ASSERT (v1.getBuffer() == 0);

        v1.set (table, ARRAY_SIZE(table));
        CPPUNIT_ASSERT (v1.size() == ARRAY_SIZE(table));
        for (size_t i=0; i<ARRAY_SIZE(table); i++)  {  CPPUNIT_ASSERT (v1[i] == table[i]);  }

        v1.set (0, 0);
        CPPUNIT_ASSERT (v1.size() == 0);
        CPPUNIT_ASSERT (v1.getBuffer() == 0);
    }

    /********************************************************************************/
    void vector_check2 ()
    {
        char table[] = { 1, 2, 3, 5, 8, 13, 21, 34, 55, 89};

        /** We create a reference vector and use it locally. */
        Vector<char>* ref = new Vector<char> (ARRAY_SIZE(table));
        ref->use ();

        /** We put some values into the vector. */
        for (size_t i=0; i<ARRAY_SIZE(table); i++)  { (*ref)[i] = table[i]; }

        /** We create other vectors referencing the previous one. */
        Vector<char> v1;  v1.setRef (ref, 0, 2);
        Vector<char> v2;  v2.setRef (ref, 2, 2);
        Vector<char> v3;  v3.setRef (ref, 4, 2);
        Vector<char> v4;  v4.setRef (ref, 6, 2);
        Vector<char> v5;  v5.setRef (ref, 8, 2);

        /** We release locally the referred vector. Since it is referred by the other vectors,
         * it won't be deleted by the release here but when the last vector will release it. */
        ref->forget ();

        /** Now we check the content of the other vectors. */
        CPPUNIT_ASSERT (v1.size() == 2);    for (size_t i=0; i<v1.size(); i++)  {  CPPUNIT_ASSERT (v1[i] == table[i+0]);  }
        CPPUNIT_ASSERT (v2.size() == 2);    for (size_t i=0; i<v2.size(); i++)  {  CPPUNIT_ASSERT (v2[i] == table[i+2]);  }
        CPPUNIT_ASSERT (v3.size() == 2);    for (size_t i=0; i<v3.size(); i++)  {  CPPUNIT_ASSERT (v3[i] == table[i+4]);  }
        CPPUNIT_ASSERT (v4.size() == 2);    for (size_t i=0; i<v4.size(); i++)  {  CPPUNIT_ASSERT (v4[i] == table[i+6]);  }
        CPPUNIT_ASSERT (v5.size() == 2);    for (size_t i=0; i<v5.size(); i++)  {  CPPUNIT_ASSERT (v5[i] == table[i+8]);  }
    }

    /********************************************************************************/

    /** \brief Check reference of a reference.
     *
     * Note: we instantiate ref1 and ref2 but don't take a local reference (with 'use')
     *  => they are referred by ref3 (directly and indirectly) and so should be deleted when needed.
     */
    void vector_check3 ()
    {
        char table[] = { 1, 2, 3, 5, 8, 13, 21, 34, 55, 89};

        /** We create a reference vector. */
        Vector<char>* ref1 = new Vector<char> (ARRAY_SIZE(table));
        for (size_t i=0; i<ARRAY_SIZE(table); i++)  { (*ref1)[i] = table[i]; }

        /** We create a reference vector. */
        Vector<char>* ref2 = new Vector<char> ();
        ref2->setRef (ref1, 3, 5);  // should hold 5, 8, 13, 21, 34
        CPPUNIT_ASSERT (ref2->size() == 5);

        /** We create a reference vector and use it locally. */
        Vector<char> ref3;
        ref3.setRef (ref2, 1, 3);  // should hold 8, 13, 21
        CPPUNIT_ASSERT (ref3.size() == 3);
        CPPUNIT_ASSERT (ref3[0] == 8);
        CPPUNIT_ASSERT (ref3[1] == 13);
        CPPUNIT_ASSERT (ref3[2] == 21);
    }

    /********************************************************************************/
    void parser_check1_aux (IOptionsParser* parser, const string& str, bool ok, size_t nbProps, const string& check)
    {
        IProperties* props = 0;
        try
        {
            props = parser->parseString (str);

            /** We can't be here. */
            CPPUNIT_ASSERT (ok);
        }
        catch (OptionFailure& e)
        {
            CPPUNIT_ASSERT (!ok);
        }

        if (ok)
        {
            CPPUNIT_ASSERT (props != 0);

            if (props->getKeys().size() > 0)
            {
                stringstream ss;
                XmlDumpPropertiesVisitor xml (ss, false, false);
                props->accept (&xml);

                CPPUNIT_ASSERT (check == ss.str());
            }

            CPPUNIT_ASSERT (props->getKeys().size() == nbProps);
        }
    }

    void parser_check1 ()
    {
        OptionsParser* parser1 = new OptionsParser ("dsk");
        parser1->push_back (new OptionOneParam ("-a",  "option a", false, "1"));

        OptionsParser* parser2 = new OptionsParser ("bloom");
        parser2->push_back (new OptionOneParam ("-b",  "option b", true));

        OptionsParser* parser3 = new OptionsParser ("debloom");
        parser3->push_back (new OptionOneParam ("-c",  "option c", false));

        OptionsParser* composite = new OptionsParser ("dbgh5");
        LOCAL (composite);

        composite->push_back (parser1);
        composite->push_back (parser2);
        composite->push_back (parser3);

        parser_check1_aux (parser1, "-a 1", true,  1, "<-a>1</-a>");
        parser_check1_aux (parser1, "-b 2", false, 1, "<-a>1</-a>");

        parser_check1_aux (parser2, "-a 3", false, 1, "");
        parser_check1_aux (parser2, "-b 4", true,  1, "<-b>4</-b>");

        parser_check1_aux (parser3, "-c 5", true, 1, "<-c>5</-c>");

        parser_check1_aux (composite, "-a 5 -b 6 -c 7", true,  3, "<-a>5</-a><-b>6</-b><-c>7</-c>");
        parser_check1_aux (composite, "-a 8 -b 9",      true,  2, "<-a>8</-a><-b>9</-b>");
        parser_check1_aux (composite, "-a 10 -c 11",    false, 0, "");
    }

    /********************************************************************************/
    void parser_check2 ()
    {
        {
            OptionsParser parser ("test1");
            parser.push_back (new OptionNoParam ("-a",  "option a", false));
            parser.push_back (new OptionNoParam ("-b",  "option b", false));
            try  {
                parser.parseString ("-a -b");
                CPPUNIT_ASSERT(parser.saw("-a"));
                CPPUNIT_ASSERT(parser.saw("-b"));
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(false);  }
        }
        {
            OptionsParser parser ("test2");
            parser.push_back (new OptionNoParam  ("-a",  "option a", false));
            parser.push_back (new OptionOneParam ("-b",  "option b", false));
            try  {
                parser.parseString ("-a -b foo");
                CPPUNIT_ASSERT(parser.saw("-a"));
                CPPUNIT_ASSERT(parser.saw("-b"));
                CPPUNIT_ASSERT(parser.getProperties()->getStr("-b") == "foo");
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(false);  }
        }
        {
            OptionsParser parser ("test3");
            parser.push_back (new OptionOneParam ("-a",  "option a", false));
            parser.push_back (new OptionNoParam  ("-b",  "option b", false));
            try  {
                parser.parseString ("-a foo -b");
                CPPUNIT_ASSERT(parser.saw("-a"));
                CPPUNIT_ASSERT(parser.saw("-b"));
                CPPUNIT_ASSERT(parser.getProperties()->getStr("-a") == "foo");
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(false);  }
        }
        {
            OptionsParser parser ("test4");
            parser.push_back (new OptionNoParam  ("-a",  "option a", false));
            parser.push_back (new OptionOneParam ("-b",  "option b", false));
            parser.push_back (new OptionNoParam  ("-c",  "option c", false));
            try  {
                parser.parseString ("-a -b foo -c");
                CPPUNIT_ASSERT(parser.saw("-a"));
                CPPUNIT_ASSERT(parser.saw("-b"));
                CPPUNIT_ASSERT(parser.getProperties()->getStr("-b") == "foo");
                CPPUNIT_ASSERT(parser.saw("-c"));
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(false);  }
        }
        {
            OptionsParser parser ("test5");
            parser.push_back (new OptionOneParam ("-b",  "option b", false));
            try  {
                parser.parseString ("-b");
                CPPUNIT_ASSERT(false);
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(true);  }
        }
        {
            OptionsParser parser ("test6");
            parser.push_back (new OptionNoParam ("-b",  "option b", false));
            try  {
                parser.parseString ("");
                CPPUNIT_ASSERT(true);
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(false);  }
        }
        {
            OptionsParser parser ("test7");
            try  {
                parser.parseString ("");
                CPPUNIT_ASSERT(true);
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(false);  }
        }
        {
            OptionsParser parser ("test7");
            parser.push_back (new OptionNoParam ("-b",  "option b", true));
            try  {
                parser.parseString ("");
                CPPUNIT_ASSERT(false);
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(true);  }
        }
        {
            OptionsParser parser ("test8");
            parser.push_back (new OptionNoParam  ("-a",  "option a", false));
            parser.push_back (new OptionOneParam ("-b",  "option b", true));
            try  {
                parser.parseString ("-a");
                CPPUNIT_ASSERT(false);
            }
            catch (OptionFailure& e)    {  CPPUNIT_ASSERT(true);  }
        }
    }

    /********************************************************************************/
    void stringline_check1 (void)
    {
        string s1 = "abcdefghijklmnopqrstuvwxyz";
        string s2 = s1 + s1;

        CPPUNIT_ASSERT (StringLine::format (s1).size() == StringLine::getDefaultWidth());
        CPPUNIT_ASSERT (StringLine::format (s2).size() == StringLine::getDefaultWidth());
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestMisc);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestMisc);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
