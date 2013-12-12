//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a kmer model with a given span size.
    Model<LargeInt<1> > model1 (27);

    // We get some information about the model.
    std::cout << "----------   MODEL   ----------"                << std::endl;
    std::cout << "span:             " << model1.getSpan()         << std::endl;
    std::cout << "kmer memory size: " << model1.getMemorySize()   << std::endl;
}
//! [snippet1]
