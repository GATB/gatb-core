//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*              Create a Collection and save it with the Storage class.         */
/*                                                                              */
/* This snippet shows how to use the Storage layer. In this example, we use the */
/* HDF5 layer, so we can check the result of the test with HDF5 tools (like     */
/* h5dump for instance).                                                        */
/*                                                                              */
/* NOTE: GATB provides some HDF5 tools (check 'bin' directory)                  */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
//! [snippet1_storage]
    // We create a Storage product "foo" in HDF5 format
    Storage* storage = StorageFactory(STORAGE_HDF5).create ("foo", true, false);
//! [snippet1_storage]

    // We use locally this object (means that it should be automatically deleted when
    // leaving the enclosing instructions block).
    LOCAL (storage);

    // Shortcut: we get the root of this Storage object
    Group& root = storage->root();

    // We get a collection of native integer from the storage.
    Collection<NativeInt64>& myIntegers = root.getCollection<NativeInt64> ("myIntegers");

    // We add some entries into the collection
    myIntegers.insert (1);
    myIntegers.insert (2);
    myIntegers.insert (3);
    myIntegers.insert (5);
    myIntegers.insert (8);

    // We flush the collection to be sure to save the content properly.
    myIntegers.flush();

    // Now, you can see the content of the collection by launching the following command:  h5dump foo.h5
    // You should get something like this:

    // HDF5 "foo.h5" {
    // GROUP "/" {
    // DATASET "myIntegers" {
    //      DATATYPE  H5T_STD_U8LE
    //      DATASPACE  SIMPLE { ( 5 ) / ( H5S_UNLIMITED ) }
    //      DATA {
    //      (0): 1, 2, 3, 5, 8
    //      }
    //  }
    // }
    // }
}
//! [snippet1]
