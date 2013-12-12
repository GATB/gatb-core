//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a kmer model with a given span size.
    Model <LargeInt<1> > model (3);

    // We compute the kmer for a given sequence
    LargeInt<1> kmer = model.codeSeed ("CAT", Data::ASCII);
    std::cout << "kmer is " << kmer << std::endl;
}
//! [snippet1]
