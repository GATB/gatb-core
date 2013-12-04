//! [snippet1]
// We include what we need for the test
#include <gatb/tools/designpattern/impl/IteratorWrappers.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::tools::dp::impl;

int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values1[] = {1,2,3,5,8,13,21,34};
    list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

    // We declare a STL list with some values.
    float values2[] = {0.5, 3.1415, 2.71};
    list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]) );

    // We declare two iterators on the two lists.
    ListIterator<int>   it1 (l1);
    ListIterator<float> it2 (l2);

    // We declare a Cartesian iterator on the two iterators.
    CartesianIterator<int,float> it (it1, it2);

    // We iterate the Cartesian product of the two lists
    for (it.first(); !it.isDone(); it.next())
    {
        cout << it->first << " -- " << it->second << endl;
    }
}
//! [snippet1]
