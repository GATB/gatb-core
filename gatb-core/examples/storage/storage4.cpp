//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                  Read two Collections from a Storage file.                   */
/*                                                                              */
/*  This snippet reads items created during 'storage2'                          */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We load a Storage product "foo" in HDF5 format
    // It should have been created with the storage2 snippet
    Storage* storage = StorageFactory(STORAGE_HDF5).load ("foo");
    LOCAL (storage);

    // Shortcut: we get the root of this Storage object
    Group& root = storage->root();

    // We get two groups from the root
    Group& group1 = root.getGroup("group1");
    Group& group2 = root.getGroup("group2");

    // We iterate the two collections with a lambda expression. Note that we use lambda expressions here.
    group1.getCollection<NativeInt64> ("integers").iterate ([] (const NativeInt64& n)  {  cout << n << endl;  });
    cout << endl;
    group2.getCollection<NativeInt64> ("integers").iterate ([] (const NativeInt64& n)  {  cout << n << endl;  });
}
//! [snippet1]
