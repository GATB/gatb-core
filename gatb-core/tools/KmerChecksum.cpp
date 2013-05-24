#include <gatb/system/impl/System.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/kmer/impl/Model.hpp>

#include <iostream>
#include <vector>

using namespace std;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::misc;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/

int main (int argc, char* argv[])
{

    if (argc < 2)
    {
        cerr << "you must provide at least 1 arguments. Arguments are:" << endl;
        cerr << "   1) uri of the file"  << endl;
        return EXIT_FAILURE;
    }

    string filename = argv[1];

    u_int64_t nbItems  = 0;
    kmer_type checksum = 0;

    IteratorFile<kmer_type> it (filename);
    for (it.first(); !it.isDone(); it.next())
    {
        checksum = checksum + *it;
        nbItems ++;
    }

    cout << "FOUND " << nbItems << " WITH CHECKSUM " << checksum << endl;

    return EXIT_SUCCESS;
}
