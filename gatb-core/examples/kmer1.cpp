//! [snippet1]
// We include what we need for the test
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/math/ttmath/ttmath.h>
#include <iostream>

// We use the required packages
using namespace std;
using namespace gatb::core::kmer::impl;

int main (int argc, char* argv[])
{
    // We declare a kmer model with a given span size.
    KmerModel model1 (27);

    // We get some information about the model.
    cout << "----------  MODEL 1  ----------"                << endl;
    cout << "span:             " << model1.getSpan()         << endl;
    cout << "kmer memory size: " << model1.getMemorySize()   << endl;

//    // We declare another kmer model with a given span size.
//    Model <ttmath::UInt<80> > model2 (80);
//
//    // We get some information about the model.
//    cout << "----------  MODEL 2  ----------"                << endl;
//    cout << "span:             " << model2.getSpan()         << endl;
//    cout << "kmer memory size: " << model2.getMemorySize()   << endl;
}
//! [snippet1]
