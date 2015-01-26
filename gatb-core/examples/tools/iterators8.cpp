//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

// Shortcut
typedef std::pair<int,float> Type;

// We define our filtering functor
struct FilterFunctor  {  bool operator ()  (const Type& val)   {  return val.first > val.second; } };

/********************************************************************************/
/*                              Mixing iterators                                */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values1[] = {1,2,3,5,8,13,21,34};
    std::list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

    // We declare a STL list with some values.
    float values2[] = {0.5, 3.1415, 12.71, -1.51, 4.11, -11.3 };
    std::list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]));

    // We declare our 'mix' iterator
    FilterIterator<Type,FilterFunctor> it (
        new PairedIterator<int,float> (
            new ListIterator<int>   (l1),
            new ListIterator<float> (l2)
        ),
        FilterFunctor()
    );

    // We iterate the pairs of the two lists
    for (it.first(); !it.isDone(); it.next())
    {
        std::cout << it->first << " -- " << it->second << std::endl;
    }
}
//! [snippet1]
