//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

/********************************************************************************/
/*            Iteration of a Cartesian product of two iterators                 */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values1[] = {1,2,3,5,8,13,21,34};
    std::list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

    // We declare a STL list with some values.
    float values2[] = {0.5, 3.1415, 2.71};
    std::list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]) );

    // We declare two iterators on the two lists.
    ListIterator<int>   it1 (l1);
    ListIterator<float> it2 (l2);

    // We declare a Cartesian iterator on the two iterators.
    ProductIterator<int,float> it (it1, it2);

    // We iterate the Cartesian product of the two lists
    for (it.first(); !it.isDone(); it.next())
    {
        std::cout << it->first << " -- " << it->second << std::endl;
    }
}
//! [snippet1]
