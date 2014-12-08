//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*         Multithreaded iteration without modification of a shared resource.   */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We get the number of cores to be used.  If we don't give any number,
    // we set to 0 which implies the usage of all available cores
    size_t nbCores = (argc >=2 ? atoi(argv[1]) : 0);

    // We create an iterator over an integer range
    int nmax = 1000;
    Range<int>::Iterator it (1,nmax);

    // We create a dispatcher configured for 'nbCores' cores.
    // The second argument tells how many consecutive values will be received by
    // each thread.
    Dispatcher dispatcher (nbCores, 1);

    // In this example, we have a different approach: we won't modify the same
    // shared integer value. Instead, each thread will use its own local integer
    // and at the end, all the local sums will be summed into the final one.
    // By doing this, we don't have any more concurrent accesses issues.

    //! [snippet5_threadobject]
    // In order to ease this approach, we use a ThreadObject object. Such an object
    // will provide local sums for each executing thread. After the iteration, it also
    // provides a mean to get all the local sums and modify the global sum accordingly.
    ThreadObject<int> sum;

    // We iterate our range.
    dispatcher.iterate (it, [&] (int i)
    {
        // We retrieve the local sum for the current executing thread with 'sum()'
        // Note that this block instruction still doesn't refer explicit thread
        // management; this is hidden through the () operator of the ThreadObject class.
        sum() += i;
    });

    // We retrieve all local sums through the 'foreach' method.
    // This loop is done 'nbCores' times.
    sum.foreach ([&] (int localSum)
    {
        // Here, the *s expression represents the global integer; we add to it the
        // current 'localSum' computed by one of the threads.
        *sum += localSum;
    });
    //! [snippet5_threadobject]

    cout << "sum=" << *sum << "  (result should be " << nmax*(nmax+1)/2 << ")" << endl;

    // In brief, the ThreadObject is a facility to avoid concurrent accesses to a shared
    // resource. It encapsulates all the tedious management of local resources and final
    // result aggregation.
}
//! [snippet1]
