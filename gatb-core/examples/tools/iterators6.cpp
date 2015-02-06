//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

/********************************************************************************/
// We a define a functor that will be called during iteration for filtering odd items.
struct FilterFunctor  {  bool operator ()  (int& val)   {  return val%2 == 0; } };

/********************************************************************************/
/*                           Iteration with filtering                           */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    std::list<int> l (values, values + sizeof(values)/sizeof(values[0]));

    // We declare a functor for filtering items.
    FilterFunctor filter;

    // We declare an iterator over list entries with filtering out odd values.
    FilterIterator<int,FilterFunctor> it (new ListIterator<int> (l), filter);

    // We iterate the truncated list
    for (it.first(); !it.isDone(); it.next())
    {
        // We should have only even values here.
        std::cout << *it << std::endl;
    }
}
//! [snippet1]
