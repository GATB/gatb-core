//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <iostream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                  Associate metadata to HDF5 collections                      */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a Storage product "foo" in HDF5 format
    Storage* storage = StorageFactory(STORAGE_HDF5).create ("foo", true, false);

    // We use locally this object (means that it should be automatically deleted when
    // leaving the enclosing instructions block).
    LOCAL (storage);

    // Shortcut: we get the root of this Storage object
    Group& root = storage->root();

    // We get a collection of native integer from the storage.
    Collection<NativeInt64>& myIntegers = root.getCollection<NativeInt64> ("myIntegers");

    // We associate a custom data to our collection
    myIntegers.addProperty ("myData", "test_%d", 147);

    // We can retrieve later this metadata with getProperty
    cout << "metadata is " << myIntegers.getProperty("myData") << endl;

    // You can dump such values with h5dump:
    //      h5dump -a myIntegers/myData foo.h5
}
//! [snippet1]
