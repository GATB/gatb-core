//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <list>
#include <iostream>

//! [snippet1_SubjectIterator]

/********************************************************************************/
// We a define a functor that will be called during iteration
struct ProgressFunctor : public IteratorListener
{
    // we receive 'current' amount of done iterations since last call
    void inc  (u_int64_t current)   {  std::cout << ".";  }
};

/********************************************************************************/
/*                     Iteration with progress information                      */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a STL list with some values.
    int values[] = {1,2,3,5,8,13,21,34};
    int valuesLen = sizeof(values)/sizeof(values[0]);
    std::list<int> l (values, values + valuesLen);

    // We create an iterator on the list.
    ListIterator<int>* itList = new ListIterator<int> (l);

    // We declare an iterator that will send progress status every 3 iterations.
    // Note that it refers a ListIterator instance given as constructor parameter.
    SubjectIterator<int> itNotif (itList, 3);

    //  We create some listener to be notified about progress and attach it to the iterator.
    itNotif.addObserver (new ProgressFunctor ());

    // We iterate the truncated list
    for (itNotif.first(); !itNotif.isDone(); itNotif.next())
    {
        // We can do something in the loop, but there is nothing here about progress management
    }
}

//! [snippet1_SubjectIterator]

//! [snippet1]
