//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*              Create 2 Collections and save them with the Storage class.      */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a Storage product "foo" in HDF5 format
    Storage* storage = StorageFactory(STORAGE_HDF5).create ("foo", true, false);
    LOCAL (storage);

    // Shortcut: we get the root of this Storage object
    Group& root = storage->root();

    // We get two groups from the root
    Group& group1 = root.getGroup("group1");
    Group& group2 = root.getGroup("group2");

    // We get two collections of native integer from the two groups.
    // Note that we can use the same 'integers' name since they will be located in
    // two different groups.
    Collection<NativeInt64>& integers1 = group1.getCollection<NativeInt64> ("integers");
    Collection<NativeInt64>& integers2 = group2.getCollection<NativeInt64> ("integers");

    // We add some entries into the collection
    integers1.insert (1);
    integers1.insert (2);
    integers1.insert (3);
    integers2.insert (5);
    integers2.insert (8);

    // We flush the collections to be sure to save the content properly.
    integers1.flush();
    integers2.flush();

    // Now, you can see the content of the collections by launching the following command:  h5dump foo.h5
    // You should get something like this:

    // HDF5 "foo.h5" {
    // GROUP "/" {
    //    GROUP "group1" {
    //       DATASET "integers" {
    //          DATATYPE  H5T_STD_U8LE
    //          DATASPACE  SIMPLE { ( 3 ) / ( H5S_UNLIMITED ) }
    //          DATA {
    //          (0): 1, 2, 3
    //          }
    //       }
    //    }
    //    GROUP "group2" {
    //       DATASET "integers" {
    //          DATATYPE  H5T_STD_U8LE
    //          DATASPACE  SIMPLE { ( 2 ) / ( H5S_UNLIMITED ) }
    //          DATA {
    //          (0): 5, 8
    //          }
    //       }
    //    }
    // }
    // }
}
//! [snippet1]
