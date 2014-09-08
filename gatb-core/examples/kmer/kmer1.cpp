//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how instantiate the Model class and how to get info       */
/* about it.                                                                    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We declare a default kmer model
    Kmer<>::ModelCanonical model;

    // We get some information about the model.
    std::cout << "----------   MODEL   ----------"               << std::endl;
    std::cout << "span:             " << model.getSpan()         << std::endl;
    std::cout << "kmer size:        " << model.getKmerSize()     << std::endl;
    std::cout << "kmer memory size: " << model.getMemorySize()   << std::endl;
}
//! [snippet1]
