//! [snippet1]
// We include what we need for the test
#include <gatb/kmer/impl/Model.hpp>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

int main (int argc, char* argv[])
{
    // We declare a kmer model with a given span size.
    Model<u_int64_t> model (3);

    // We compute the kmer for a given sequence
    u_int64_t kmer = model.codeSeed ("CAT", KMER_DIRECT);
    cout << "kmer is " << kmer << endl;

    // We compute the next kmer on the right by adding one nucleotide.
    kmer = model.codeSeedRight (kmer, 'T', KMER_DIRECT);
    cout << "new kmer is " << kmer << endl;
}
//! [snippet1]
