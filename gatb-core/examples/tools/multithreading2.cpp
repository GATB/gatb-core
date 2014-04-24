//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*         Multithreaded iteration and modification of a shared resource.       */
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
    // each thread. The second argument tells how to group items per thread (set
    // here to 1 to emphasize concurrent access issue).
    Dispatcher dispatcher (nbCores, 1);

    // The idea here is to sum the integers of our range with an iteration.
    // (Note: we know that the result is N*(N+1)/2)
    int sum1=0, sum2=0;

    //////////////////////////////////////////////////
    // First iteration: WRONG WAY
    //////////////////////////////////////////////////
    // Our first attempt is to use an integer variable to sum the iterated value.
    // This variable will be shared by all the threads and, since they access to it
    // without caution wrt concurrent accesses, the sum result should be wrong (unless
    // you use one core only)
    dispatcher.iterate (it, [&] (int i)  {  sum1 += i;  });

    //////////////////////////////////////////////////
    // Second iteration: CORRECT WAY
    //////////////////////////////////////////////////
    // As previously, our second attempt will share the same integer variable.
    // But now, we take care about concurrent accesses with the use of the
    // __sync_fetch_and_add intrinsic instruction. This instruction ensures that
    // the shared integer can be modified by only one thread at one time.
    dispatcher.iterate (it, [&] (int i)  {  __sync_fetch_and_add (&sum2, i);  });

    //////////////////////////////////////////////////
    // CONCLUSION
    //////////////////////////////////////////////////
    cout << "First iteration:  sum=" << sum1 << "  (result should be " << nmax*(nmax+1)/2 << ")" << endl;
    cout << "Second iteration: sum=" << sum2 << "  (result should be " << nmax*(nmax+1)/2 << ")" << endl;

    // Parallelization of Iterator is pretty simple with the Dispatcher class.
    // Moreover, usage of lambda expressions make the whole thing easy to write.
    // Note that the instruction block of the lambda expression doesn't even know that
    // it may be executed in different threads. In other words, the block doesn't refer
    // any stuff related to thread management; it just receives one of the item of the
    // iteration and process some action on it.

    // IMPORTANT ! As we have seen here, the user has to be aware that a shared resource (one
    // integer here) can be modified by several threads at the same time, so the user must use
    // some kind of synchronization for modifying the shared resource. We will see in other
    // examples that GATB provides mechanisms for this purpose.
}
//! [snippet1]
