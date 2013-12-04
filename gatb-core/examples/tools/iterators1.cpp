//! [snippet1]
// We include what we need for the test
#include <gatb/tools/designpattern/impl/IteratorWrappers.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::tools::dp::impl;

int main (int argc, char* argv[])
{
    size_t nbItems = 0;

    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    list<int> l (values, values + sizeof(values)/sizeof(values[0]) );

    // We create an iterator on this STL list
    ListIterator<int> it (l);

    // We iterate the list through the iterator.
    for (it.first(); !it.isDone(); it.next())
    {
        cout << *it << endl;
    }
}
//! [snippet1]
