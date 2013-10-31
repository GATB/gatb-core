//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <string.h>

// We use the required packages
using namespace std;

int main (int argc, char* argv[])
{
    const char* seq = "CATTGATAGTGG";

    // We declare a kmer model with a given span size.
    Model<LargeInt<1> > model (3);

    // We declare an iterator on a given sequence.
    Model<LargeInt<1> >::Iterator it (model);

    // We set the data from which we want to extract kmers.
    Data data (Data::ASCII);
    data.set ((char*)seq, strlen(seq));

    it.setData (data);

    // We iterate the kmers.
    for (it.first(); !it.isDone(); it.next())
    {
        cout << *it << endl;
    }
}
//! [snippet1]
