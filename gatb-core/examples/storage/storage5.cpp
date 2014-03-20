//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <iostream>
#include <memory>

// We use the required packages
using namespace std;

/********************************************************************************/
/*                  Iterate solid kmers from a HDF5 file.                       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We check that the user provides a graph URL (supposed to be in HDF5 format).
    if (argc < 2)
    {
        std::cerr << "You must provide a HDF5 file." << std::endl;
        return EXIT_FAILURE;
    }

    // We get a handle on the HDF5 storage object.
    // Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
    auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (argv[1]));

    // We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
    Collection<Kmer<>::Count>& solidKmers = storage->getGroup("dsk").getCollection<Kmer<>::Count> ("solid");

    // Note the type of the collection: Kmer<span>::Count
    // If we don't specify 'span', we get kmers of size 32

    // We iterate (through a lambda expression) the solid kmers from the retrieved collection
    size_t nbKmers = 0;
    solidKmers.iterate ([&] (const Kmer<>::Count& count)
    {
        // We dump the solid kmer value and its abundance
        cout << "[" << ++nbKmers << "]  " << count.value << "  "  << count.abundance << endl;
    });

    // We can also retrieve information (as an XML string) about the construction of the solid kmers
    cout << solidKmers.getProperty("properties") << endl;

    // We can access each of these information through a Properties object
    Properties props;
    props.readXML (solidKmers.getProperty("properties"));

    // Now, we can for instance get the kmer size (as an integer)
    cout << "kmer size:      " << props.getInt ("kmer_size")      << endl;
    cout << "nb solid kmers: " << props.getInt ("kmers_nb_solid") << endl;
}
//! [snippet1]
