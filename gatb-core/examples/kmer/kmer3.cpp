//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
/*                              Kmer management                                 */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We configure a data object with a sequence (in ASCII format)
    Data data ((char*)"CATTGATAGTGG");

    // We declare a kmer model with a given span size.
    Model<LargeInt<1> > model (3);

    // We declare an iterator on a given sequence.
    Model<LargeInt<1> >::Iterator it (model);

    // We configure the iterator with our sequence
    it.setData (data);

    // We iterate the kmers.
    for (it.first(); !it.isDone(); it.next())
    {
        std::cout << *it << std::endl;
    }
}
//! [snippet1]
