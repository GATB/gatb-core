//! [snippet1]
// We include what we need for the test
#include <gatb/tools/designpattern/impl/IteratorWrappers.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

// We a define a functor that will be called during iteration for filtering odd items.
struct FilterFunctor  {  bool operator ()  (int& val)   {  return val%2 == 0; } };

int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    int valuesLen = sizeof(values)/sizeof(values[0]);
    list<int> l (values, values + valuesLen);

    // We declare an iterator on the list.
    ListIterator<int> listIt (l);

    // We declare a functor for filtering items.
    FilterFunctor filter;

    // We declare an iterator over list entries with a filter on odd values.
    FilterIterator<int,FilterFunctor> it (listIt, filter);

    // We iterate the truncated list
    for (it.first(); !it.isDone(); it.next())
    {
        // We should have only even values here.
        cout << *it << endl;
    }
}
//! [snippet1]
