//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

/********************************************************************************/
/*            Iteration of a two iterators by pairs of items                    */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values1[] = {13,5,34};
    std::list<int> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

    // We declare a STL list with some values.
    float values2[] = {0.5, 3.1415, 2.71};
    std::list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]) );

    // We declare a paired iterator on the two iterators.
    PairedIterator<int,float> it (new ListIterator<int>(l1), new ListIterator<float>(l2));

    // We iterate the pairs of the two lists
    for (it.first(); !it.isDone(); it.next())
    {
        std::cout << it->first << " -- " << it->second << std::endl;
    }
}
//! [snippet1]
