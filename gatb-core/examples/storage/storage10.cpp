// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <iostream>

// We use the required packages
using namespace std;

/********************************************************************************/
/*  This snippet shows how to use binary input/output streams with files.       */
/*  so it's the same as storage8 but with files instead of hdf5.                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    float table[] = { 0.577, 3.1415, 2.71 };

    Storage* storage = StorageFactory(STORAGE_FILE).create ("foo", true, false);

    // We use locally this object (means that it should be automatically deleted when
    // leaving the enclosing instructions block).
    LOCAL (storage);

    // Shortcut: we get the root of this Storage object
    Group& root = storage->root();

    // We get an output stream in a C++ style
    Storage::ostream os (root, "data");

    // We write some information in this stream
    os.write (reinterpret_cast<char const*>(table), sizeof(table));

    // We have to flush the stream in order to be sure everything is ok
    os.flush();

    // We get a handle on the collection where we put our data
    // Note: the collection is typed as NativeInt8, meaning we get binary data
    Collection<NativeInt8>& dataCollection = root.getCollection<NativeInt8> ("data");

    // We get the number of items in the collection.
    size_t nbItems = dataCollection.getNbItems() / sizeof(table[0]);
    cout << "nb items : " << nbItems << endl;

    // Now we declare an input stream on the collection
    Storage::istream is (root, "data");

    // We want to read the data, we first need to have a buffer for this
    float* buffer = new float [nbItems];

    // We read the data from the input stream
    is.read (reinterpret_cast<char*>(buffer), nbItems*sizeof(float));

    // We check that we read correct values.
    cout << "check : " << (memcmp(buffer, table, nbItems*sizeof(float)) == 0) << endl;

    // cleanup
    delete[] buffer;
}

