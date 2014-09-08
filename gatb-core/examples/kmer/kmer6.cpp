//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* A Model is defined:                                                          */
/*      - the maximum kmer size supported (the template argument)               */
/*      - the actual kmer size used (the constructor argument)                  */
/*                                                                              */
/* This snippet shows that the actual size can't exceed the max size.           */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a kmer model that supports size up to 32
    // Unfortunately, we try to have a model with kmer of size 51
    // => we should get an exception
    try
    {
        Kmer<32>::ModelCanonical model (51);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    try
    {
        Kmer<64>::ModelCanonical model (51);
        std::cout << "It's OK..." << std::endl;
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
