//! [snippet1]
// GATB-Core Online Tutorial 

/********************************************************************************/
/*                                                                              */
/*  Multithreading: iteration and modification of a shared resource.            */
/*                                                                              */
/*  We illustrate here how to write to a file, which is definitively not an     */ 
/*  operation that be done using concurrent accesses.                           */
/*                                                                              */
/********************************************************************************/
#include <gatb/gatb_core.hpp>

#include <fstream>
using namespace std;

// We define a functor that will be automatically cloned by the dispatcher.
// nb clones = nb cores!
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

  // We create the synchronizer that will be shared between all concurrent
  // threads that are going to be create by the Dispatcher, just below.
  ISynchronizer* synchro = System::thread().newSynchronizer();

  // We create a dispatcher configured for 'nbCores' cores.
  Dispatcher dispatcher (nbCores, 1);

  // We iterate over the range of integers.  
  dispatcher.iterate (it, Functor(synchro,file));

  // We close the file
  file.close();

  // We get rid of the synchronizer
  delete synchro;
}
//! [snippet1]
