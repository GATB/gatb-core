//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

/********************************************************************************/
/*                     Iteration of a part of an iterator                       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    int valuesLen = sizeof(values)/sizeof(values[0]);
    std::list<int> l (values, values + valuesLen);

    // We declare one iterator on the list
    ListIterator<int> it (l);

    // We declare a truncated iterator for the list iterator.
    TruncateIterator<int> itTrunc (it, valuesLen/2);

    // We iterate the truncated list
    for (itTrunc.first(); !itTrunc.isDone(); itTrunc.next())
    {
        std::cout << *itTrunc << std::endl;
    }
}
//! [snippet1]
