//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <vector>
#include <iostream>

/********************************************************************************/
/*            Iteration of a two iterators by pairs of items                    */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    float values1[] = {13,5,34};
    std::list<float> l1 (values1, values1 + sizeof(values1)/sizeof(values1[0]) );

    // We declare a STL list with some values.
    float values2[] = {0.5, 3.1415, 2.71};
    std::list<float> l2 (values2, values2 + sizeof(values2)/sizeof(values2[0]) );

    std::vector<Iterator<float>*> iterators;
    iterators.push_back (new ListIterator<float>(l1));
    iterators.push_back (new ListIterator<float>(l2));

    // We declare a composite iterator for the two iterators.
    CompositeIterator<float> it (iterators);

    // We iterate the pairs of the two lists
    for (it.first(); !it.isDone(); it.next())
    {
        std::cout << it.item() << std::endl;
    }
}
//! [snippet1]
