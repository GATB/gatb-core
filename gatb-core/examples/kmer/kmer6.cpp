//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

int main (int argc, char* argv[])
{
    // We declare a kmer model that supports size up to 32
    // Unfortunately, we try to have a model with kmer of size 51
    // => we should get an exception
    try
    {
        Kmer<32>::Model model (51);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
