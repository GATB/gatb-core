//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how to get kmer information from a raw ascii sequence.    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We define some nucleotides sequence.
    const char* seq = "CTACGAATA";

    // We declare a kmer model with kmer size big enough to represent our sequence.
    Kmer<>::Model model (strlen(seq));

    // We compute the kmer for a given sequence
    Kmer<>::Type kmer = model.codeSeed (seq, Data::ASCII);

    std::cout << "kmer  value is: " << kmer                 << std::endl;
    std::cout << "kmer string is: " << model.toString(kmer) << std::endl;
}
//! [snippet1]
