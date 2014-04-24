//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <fstream>
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

    // We open a file. This will be our shared resource between threads.
    fstream file ("out", std::fstream::out);

    // For our file, we can't use intrinsics like we did for integer addition,
    // so we need a general synchronization mechanism that will be shared by the threads.
    ISynchronizer* synchro = System::thread().newSynchronizer();

    // We create a dispatcher configured for 'nbCores' cores.
    Dispatcher dispatcher (nbCores, 1);

    // We iterate the range
    dispatcher.iterate (it, [&] (int i)
    {
        // We use a helper class that will protect the full containing instruction block
        // against concurrent access. Note it uses our shared synchro object.
        // We don't have to do the tandem lock/unlock, a single LocalSynchronizer
        // declaration will protect the containing block. This may be useful because
        // if the user forget to call the 'unlock' method, it would block the full
        // program execution for ever.
        LocalSynchronizer sync (synchro);

        // We dump the current integer into the file
        file << i << endl;
    });

    // We close the file
    file.close();

    // We get rid of the synchronizer
    delete synchro;
}
//! [snippet1]
