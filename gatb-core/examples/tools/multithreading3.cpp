//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <fstream>
using namespace std;

/********************************************************************************/
/*         Multithreaded iteration and modification of a shared resource.       */
/********************************************************************************/

// We define a functor that will be cloned by the dispatcher
struct Functor
{
    ISynchronizer* synchro;    fstream& file;
    Functor (ISynchronizer* synchro, fstream& file)  : synchro(synchro), file(file) {}

    void operator() (int i)
    {
        // We lock the synchronizer
        synchro->lock ();

        // We dump the current integer into the file
        file << i << endl;

        // We unlock the synchronizer
        synchro->unlock ();
    }
};

/********************************************************************************/
int main (int argc, char* argv[])
{
    // We get the number of cores to be used.  If we don't give any number,
    // we set to 0 which implies the usage of all available cores
    size_t nbCores = (argc >=2 ? atoi(argv[1]) : 0);

    // We create an iterator over an integer range
    int nmax = 1000;
    Range<int>::Iterator it (1,nmax);

    // We open a file. This will be our shared resource between threads.
    fstream file ("out", std::fstream::out);

    // For our file, we can't use intrinsics like we did for integer addition,
    // so we need a general synchronization mechanism that will be shared by the threads.
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // We create a dispatcher configured for 'nbCores' cores.
    Dispatcher dispatcher (nbCores, 1);

    // We iterate the range.  NOTE: we could also use lambda expression (easing the code readability)
   dispatcher.iterate (it, Functor(synchro,file));

    // We close the file
    file.close();

    // We get rid of the synchronizer
    delete synchro;
}
//! [snippet1]
