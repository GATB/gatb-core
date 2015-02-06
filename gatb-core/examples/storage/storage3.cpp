//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                  Read a Collection from a Storage file.                      */
/*                                                                              */
/*  This snippet reads items created during 'storage1'                          */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We load a Storage product "foo" in HDF5 format
    // It must have been created with the storage1 snippet
//! [snippet1_storage]
    Storage* storage = StorageFactory(STORAGE_HDF5).load ("foo");
//! [snippet1_storage]
    LOCAL (storage);

    // Shortcut: we get the root of this Storage object
    Group& root = storage->root();

    // We get a collection of native integer from the storage.
    Collection<NativeInt64>& myIntegers = root.getCollection<NativeInt64> ("myIntegers");

    // We create an iterator for our collection.
    Iterator<NativeInt64>* iter = myIntegers.iterator();
    LOCAL (iter);

    // Now we can iterate the collection through this iterator.
    for (iter->first(); !iter->isDone(); iter->next())  {  cout << iter->item() << endl;  }
}
//! [snippet1]
