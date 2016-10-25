//! [snippet1]

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how to iterate the kmers of a given string.               */
/*                                                                              */
/********************************************************************************/

// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>

/********************************************************************************/
//START application
int main (int argc, char* argv[])
{
  // We define a sequence of nucleotides
  const char* seq = "CATTGATAGTGGATGGT";
  std::cout << "Initial sequence: " << seq << std::endl;

  // We declare a kmer model with kmer of size 5.
  // Note that we want "direct" kmers, not the min(forward,revcomp) 
  // default behavior.
  Kmer<>::ModelDirect model (5);
  // We declare an iterator relying on that k-mer model.
  Kmer<>::ModelDirect::Iterator it (model);

  // We create a data object representation of the sequence (in ASCII format)
  Data data ((char*)seq);

  // We configure the iterator with our sequence
  it.setData (data);

  // And now, we can iterate over the kmers.
  for (it.first(); !it.isDone(); it.next())
  {
    std::cout << "kmer " << model.toString(it->value()) << ",  value " << it->value() << std::endl;
  }
}
//! [snippet1]
