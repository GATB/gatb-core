//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how to iterate the kmers of a given string.               */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We define a sequence of nucleotides
    const char* seq = "CATTGATAGTGGATGGT";
    std::cout << "Initial sequence: " << seq << std::endl;

    // We configure a data object with a sequence (in ASCII format)
    Data data ((char*)seq);

    // We declare a kmer model with kmer of size 5.
    // Note that we want "direct" kmers, not the min(forward,revcomp) default behavior.
    Kmer<>::ModelDirect model (5);

    // We declare an iterator on a given sequence.
    Kmer<>::ModelDirect::Iterator it (model);

    // We configure the iterator with our sequence
    it.setData (data);

    // We iterate the kmers.
    for (it.first(); !it.isDone(); it.next())
    {
        std::cout << "kmer " << model.toString(it->value()) << ",  value " << it->value() << std::endl;
    }
}
//! [snippet1]
