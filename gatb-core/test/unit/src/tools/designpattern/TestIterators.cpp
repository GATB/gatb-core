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

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/math/Integer.hpp>

#include <vector>
#include <string>

#include <typeinfo>

#include <boost/variant.hpp>
using namespace std;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::math;

/********************************************************************************/
namespace gatb  {  namespace tests  {
/********************************************************************************/

/** \brief Test class for miscellaneous operations
 */
class TestIterators : public Test
{
    /********************************************************************************/
    CPPUNIT_TEST_SUITE_GATB (TestIterators);

        CPPUNIT_TEST_GATB (iterators_checkNullIterator);
        CPPUNIT_TEST_GATB (iterators_checkList);
        CPPUNIT_TEST_GATB (iterators_checkCartesianIterator);
        CPPUNIT_TEST_GATB (iterators_checkCompoundIterator);
        CPPUNIT_TEST_GATB (iterators_checkTruncateIterator);
        CPPUNIT_TEST_GATB (iterators_checkCancellableIterator);
        CPPUNIT_TEST_GATB (iterators_checkPairedIterator);
        CPPUNIT_TEST_GATB (iterators_checkVariant1);
        CPPUNIT_TEST_GATB (iterators_checkVariant2);
        CPPUNIT_TEST_GATB (iterators_adaptator);

    CPPUNIT_TEST_SUITE_GATB_END();

public:
    /********************************************************************************/
    void setUp    ()  {}
    void tearDown ()  {}

    /********************************************************************************/
    /** \brief check the NULL iterator
     *
     * Test of \ref gatb::core::tools::dp::impl::NullIterator \n
     */
    void iterators_checkNullIterator ()
    {
        size_t nbItems = 0;

        /** We declare a null iterator. */
        NullIterator<int> it;
        for (it.first(); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (false);
        }
        CPPUNIT_ASSERT (nbItems == 0);
    }

    /********************************************************************************/
    /** \brief check the std::list iterator wrapper.
     *
     * Test of \ref gatb::core::tools::dp::impl::ListIterator \n
     */
    void iterators_checkList ()
    {
        size_t nbItems = 0;

        /** We declare a STL list with some values. */
        int values[] = {1,2,3,5,8,13,21,34};
        list<int> l (values, values + sizeof(values)/sizeof(values[0]) );

        /** We iterate the list through a ListIterator instance. */
        ListIterator<int> it (l);
        for (it.first(); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (*it == values[nbItems++]);
        }
        CPPUNIT_ASSERT (nbItems == sizeof(values)/sizeof(values[0]));
    }

    /********************************************************************************/
    /** \brief check the CartesianIterator
     *
     * Test of \ref gatb::core::tools::dp::impl::CartesianIterator\n
     */
    void iterators_checkCartesianIterator ()
    {
        size_t nbItems = 0;

        /** We declare a STL list with some values. */
        int values1[] = {1,2,3,5,8,13,21,34};
        list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

        /** We declare a STL list with some values. */
        float values2[] = {0.5, 3.1415, 2.71};
        list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]) );

        /** We build the pairs list that should be iterated by the Cartesian iterator. */
        vector <pair<int,float> > pairsCheck;
        for (size_t i=0; i<sizeof(values1)/sizeof(values1[0]); i++)
        {
            for (size_t j=0; j<sizeof(values2)/sizeof(values2[0]); j++)
            {
                pairsCheck.push_back (pair<int,float>(values1[i],values2[j]));
            }
        }

        /** We declare two iterators on the two lists. */
        ListIterator<int>   it1 (l1);
        ListIterator<float> it2 (l2);

        /** We declare a Cartesian iterator on the two iterators. */
        ProductIterator<int,float> it (it1, it2);
        for (it.first(); !it.isDone(); it.next())
        {
            CPPUNIT_ASSERT (it->first  == pairsCheck[nbItems].first);
            CPPUNIT_ASSERT (it->second == pairsCheck[nbItems].second);
            nbItems++;
        }
        CPPUNIT_ASSERT (nbItems == sizeof(values1)/sizeof(values1[0]) * sizeof(values2)/sizeof(values2[0]));
    }

    /********************************************************************************/
    /** We define an iterator that loops every ith character of a given string */
    class MyIterator : public Iterator<char>
    {
    public:
        MyIterator (const string& s) : _str(s), _mod(1), _idx(0) {}
        void update (int modulo)  { _mod = modulo; }
        void first  () { _idx  = _mod; }
        void next   () { _idx += _mod; }
        bool isDone () { return _idx >= (int)_str.size(); }
        char& item  () { return _str[_idx]; }
    private:
        string _str;
        int    _mod;
        int    _idx;
    };

    /** We need a functor for updating the inner loop iterator. */
    struct Update { void operator() (Iterator<char>* it2, int* val)  {  static_cast<MyIterator*> (it2)->update (*val);  }};

    /** \brief check the CompoundIterator
     *
     * Test of \ref gatb::core::tools::dp::impl::CompoundIterator\n
     */
    void iterators_checkCompoundIterator ()
    {
        size_t nbItems = 0;
        const char* str = "abcdefghijklmnopqrstuvwxyz";

        /** We declare a STL list with some values. */
        int values1[] = {2, 3, 5};
        list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );
        ListIterator<int> it1 (l1);

        /** We build the check table holding the items we are supposed to iterate. */
        char checkTable[] = {
            'c','e','g','i','k','m','o','q','s','u','w','y',    // loop modulo 2
            'd','g','j','m','p','s','v','y',                    // loop modulo 3
            'f','k','p','u','z'                                 // loop modulo 5
        };

        /** We declare an iterator on our custom iterator. */
        MyIterator it2 (str);

        /** We declare a Cartesian iterator on the two iterators. */
        CompoundIterator<int,char,Update> it (it1, it2, Update());
        for (it.first(); !it.isDone(); it.next(), nbItems++)
        {
            if (nbItems < sizeof(checkTable)/sizeof(checkTable[0]))
            {
                CPPUNIT_ASSERT (*it == checkTable[nbItems]);
            }
        }

        /** We  check that we got the correct number of iterated items. */
        CPPUNIT_ASSERT (nbItems == sizeof(checkTable)/sizeof(checkTable[0]));
    }

    /********************************************************************************/
    /** \brief check the truncation of the given iterator
     *
     * Test of \ref gatb::core::tools::dp::impl::TruncateIterator \n
     */
    void iterators_checkTruncateIterator ()
    {
        int nbItems = 0;

        /** We declare a STL list with some values. */
        int values[] = {1,2,3,5,8,13,21,34};
        int valuesLen = sizeof(values)/sizeof(values[0]);
        list<int> l (values, values + valuesLen);

        /** We declare an iterator on this list and loop the items. */
        ListIterator<int> itRef (l);
        nbItems = 0;
        for (itRef.first(); !itRef.isDone(); itRef.next())  {  CPPUNIT_ASSERT (*itRef == values[nbItems++]);  }
        CPPUNIT_ASSERT (nbItems == valuesLen);

        /** We declare a truncated iterator for the list iterator. */
        TruncateIterator<int> itTrunc (itRef, valuesLen/2);
        nbItems = 0;
        for (itTrunc.first(); !itTrunc.isDone(); itTrunc.next())  {  CPPUNIT_ASSERT (*itTrunc == values[nbItems++]);  }
        CPPUNIT_ASSERT (nbItems == valuesLen/2);

        /** We declare a truncated iterator for the list iterator. Note that we make it longer than the referred one. */
        TruncateIterator<int> itTrunc2 (itRef, 2*valuesLen);
        nbItems = 0;
        for (itTrunc2.first(); !itTrunc2.isDone(); itTrunc2.next())  {  CPPUNIT_ASSERT (*itTrunc2 == values[nbItems++]);  }
        CPPUNIT_ASSERT (nbItems == valuesLen);
    }

    /********************************************************************************/
    /** \brief check the cancellation of the given iterator
     *
     * Test of \ref gatb::core::tools::dp::impl::CancellableIterator \n
     */
    void iterators_checkCancellableIterator ()
    {
        int nbItems = 0;

        /** We declare a STL list with some values. */
        int values[] = {1,2,3,5,8,13,21,34};
        int valuesLen = sizeof(values)/sizeof(values[0]);
        list<int> l (values, values + valuesLen);

        /** We declare an iterator on this list and loop the items. */
        ListIterator<int> itRef (l);
        nbItems = 0;
        for (itRef.first(); !itRef.isDone(); itRef.next())  {  CPPUNIT_ASSERT (*itRef == values[nbItems++]);  }
        CPPUNIT_ASSERT (nbItems == valuesLen);

        /** We declare a truncated iterator for the list iterator. */
        CancellableIterator<int> itCanc (itRef);
        nbItems = 0;
        for (itCanc.first(); !itCanc.isDone(); itCanc.next())  {  
            CPPUNIT_ASSERT (*itCanc == values[nbItems++]);
            if (nbItems == valuesLen/2) itCanc._cancel = true;
        }
        CPPUNIT_ASSERT (nbItems == valuesLen/2 );
    }


    /********************************************************************************/
    /** \brief check the PairedIterator
     *
     * Test of \ref gatb::core::tools::dp::impl::PairedIterator\n
     */
    void iterators_checkPairedIterator ()
    {
        size_t nbItems = 0;
        size_t i=0;

        /** We declare a STL list with some values. */
        int values1[] = {1,2,3,5,8,13,21,34};
        list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

        /** We declare a STL list with some values. */
        float values2[] = {0.5, 3.1415, 2.71, 5.87451};
        list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]) );

        /** We declare a paired iterator on the two iterators. */
        PairedIterator<int,float> it (new ListIterator<int> (l1), new ListIterator<float> (l2));
        for (it.first(); !it.isDone(); it.next(), i++)
        {
            CPPUNIT_ASSERT (it->first  == values1[i]);
            CPPUNIT_ASSERT (it->second == values2[i]);

            nbItems++;
        }

        CPPUNIT_ASSERT (nbItems > 0);
        CPPUNIT_ASSERT (nbItems == min (sizeof(values1)/sizeof(values1[0]), sizeof(values2)/sizeof(values2[0])));
    }

    /********************************************************************************/
    void iterators_checkVariant1 ()
    {
        vector<int>    v1;  v1.push_back(1);                            VectorIterator<int>    it1 (v1);
        vector<float>  v2;  v2.push_back(3.14);   v2.push_back(0.577);  VectorIterator<float>  it2 (v2);
        vector<string> v3;  v3.push_back("foo");  v3.push_back("bar");  VectorIterator<string> it3 (v3);

        IteratorVariant<VectorIterator,int,float,string> it;

        it = it3;

        size_t i=0;
        for (it.first(); !it.isDone(); it.next(), i++)
        {
        }
    }

    /********************************************************************************/
    void iterators_checkVariant2 ()
    {
		// IteratorVariant<VectorIterator,INTEGER_TYPES> it;
		//
		// vector<LargeInt<1> > v1;  v1.push_back(1);  VectorIterator<LargeInt<1> >  it1 (v1);
		// vector<LargeInt<2> > v2;  v2.push_back(2);  v2.push_back(3);  VectorIterator<LargeInt<2> > it2 (v2);

        //it = it1;
        //size_t i=0;
        //for (it.first(); !it.isDone(); it.next(), i++)  {  }
    }

    /********************************************************************************/
    struct Entry { int n; float x; };


    struct Adaptator  {  float& operator() (Entry& e)  { return e.x; }  };

    void iterators_adaptator ()
    {
        size_t i=0;

        Entry table[] = { {3,0.4}, {17,2.7}, {-2,3.14} };
        list<Entry> l (table, table + sizeof(table)/sizeof(table[0]) );

        Iterator<Entry>* it = new ListIterator<Entry> (l);
        LOCAL (it);

        i=0;
        for (it->first(); !it->isDone(); it->next(), i++)
        {
            CPPUNIT_ASSERT (it->item().n == table[i].n);
            CPPUNIT_ASSERT (it->item().x == table[i].x);
        }

        IteratorAdaptor<Entry,float,Adaptator> itAdapt (it);

        i=0;
        for (itAdapt.first(); !itAdapt.isDone(); itAdapt.next(), i++)
        {
            CPPUNIT_ASSERT (itAdapt.item() == table[i].x);
        }
    }
};

/********************************************************************************/

CPPUNIT_TEST_SUITE_REGISTRATION      (TestIterators);
CPPUNIT_TEST_SUITE_REGISTRATION_GATB (TestIterators);

/********************************************************************************/
} } /* end of namespaces. */
/********************************************************************************/
