//! [snippet1]
// We include what we need for the test
#include <gatb/kmer/impl/Model.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;
using namespace gatb::core::tools::misc;

int main (int argc, char* argv[])
{
    // We declare a kmer model with a given span size.
    KmerModel model (3);

    // We compute the kmer for a given sequence
    kmer_type kmer = model.codeSeed ("CAT", Data::ASCII);
    cout << "kmer is " << kmer << endl;
}
//! [snippet1]
