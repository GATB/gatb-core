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
    int values[] = {1,2,3,5,8,13,21,34};
    int valuesLen = sizeof(values)/sizeof(values[0]);
    list<int> l (values, values + valuesLen);

    // We declare one iterator on the list
    ListIterator<int> it (l);

    // We declare a truncated iterator for the list iterator.
    TruncateIterator<int> itTrunc (it, valuesLen/2);

    // We iterate the truncated list
    for (itTrunc.first(); !itTrunc.isDone(); itTrunc.next())
    {
        cout << *itTrunc << endl;
    }
}
//! [snippet1]
