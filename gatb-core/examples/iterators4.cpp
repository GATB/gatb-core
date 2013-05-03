//! [snippet1]
// We include what we need for the test
#include <gatb/tools/designpattern/impl/IteratorWrappers.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

//! [snippet1_SubjectIterator]

// We a define a functor that will be called during iteration
struct ProgressFunctor : public IteratorListener
{
    // we receive 'current' amount of done iterations since last call
    void inc  (u_int64_t current)   {  cout << ".";  }
};

int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    int valuesLen = sizeof(values)/sizeof(values[0]);
    list<int> l (values, values + valuesLen);

    // We create an iterator on the list.
    ListIterator<int> itList (l);

    // We declare an iterator that will send progress status every 3 iterations.
    // Note that it refers a ListIterator instance given as constructor parameter.
    SubjectIterator<int> itNotif (itList, 3);

    //  We create some listener to be notified about progress and attach it to the iterator.
    ProgressFunctor fct;    itNotif.addObserver (fct);

    // We iterate the truncated list
    for (itNotif.first(); !itNotif.isDone(); itNotif.next())
    {
        // We can do something in the loop, but there is nothing here about progress management
    }
}

//! [snippet1_SubjectIterator]

//! [snippet1]
