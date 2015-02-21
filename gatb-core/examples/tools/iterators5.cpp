//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

//! [snippet1_SubjectIterator]

/********************************************************************************/
/*                     Iteration with progress information                      */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    int valuesLen = sizeof(values)/sizeof(values[0]);
    std::list<int> l (values, values + valuesLen);

    // We declare an iterator that will send default progress status.
    // Note that we define the 'actual' iterator on the fly as first parameter of ProgressIterator
    ProgressIterator<int> it (new ListIterator<int> (l), "Iteration running", valuesLen);

    // We iterate the list
    for (it.first(); !it.isDone(); it.next())
    {
        // We force a small wait
        sleep (1);
    }
}

//! [snippet1_SubjectIterator]

//! [snippet1]
