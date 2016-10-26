//! [snippet1]
// GATB-Core Online Tutorial 

/********************************************************************************/
/*                     Multithreaded iteration of a bank.                       */
/*                                                                              */
/* Notice: usage of Dispatcher has sense only if the iterated items             */
/*         can be processed independently from each other.                      */
/*                                                                              */
/*         The point to understand with the Dispatcher is that it can           */
/*         iterate any instance of Iterator class. If you have any set of items */
/*         that can be enumerated through an Iterator implementation, then you  */
/*         can parallelize the iteration with a Dispatcher instance.            */
/*                                                                              */
/* WARNING! This snippet also shows how to use lambda expressions, so you need  */
/* to use a compiler that supports this feature.                                */
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
// START Application
int main (int argc, char* argv[])
{
  // We get the file name from the user arguments
  // Online GATB-Tutorial: this argument is automatically filled in with an 
  // appropriate file.
  if (argc < 2)
  {
  cerr << "Please, provide a sequence file." << endl;
  return EXIT_FAILURE;
  }

  // We get a handle on a bank
  BankFasta bank (argv[1]);

  // We get the number of cores to be used.  If we don't give any number,
  // we set to 0 which implies the usage of all available cores.
  size_t nbCores = (argc >=3 ? atoi(argv[2]) : 0);

  // We create a dispatcher.
  // Exercise: try to execute this code with:
  //   a. nbCores=1
  //   b. nbCores=4
  // ... then have a look at the running time
  Dispatcher dispatcher (nbCores);

  // We will count nucleotides occurrences.
  ThreadObject<int> sumA, sumC, sumG, sumT, sumN;

  // We iterate over the bank sequences. 
  // Note how we provide a bank iterator to the dispatcher.
  IDispatcher::Status status = dispatcher.iterate (bank.iterator(), [&] (const Sequence& seq)
  {
    // We use shortcut references for the different local sums. 
    // It avoids to retrieve them each time a nucleotide of the 
    // sequence is handled (see for loop below) and may give much 
    // better performance.
    int& localA = sumA();
    int& localC = sumC();
    int& localG = sumG();
    int& localT = sumT();
    int& localN = sumN();

    // We loop over the nucleotides of the current sequence.
    for (size_t i=0; i<seq.getDataSize(); i++)
    {
      switch (seq.getDataBuffer()[i])
      {
        case 'A':  localA++;  break;
        case 'C':  localC++;  break;
        case 'G':  localG++;  break;
        case 'T':  localT++;  break;
        case 'N':  localN++;  break;
      }
    }
  });

  // Finalize composition computing
  sumA.foreach ([&] (int n) { *sumA += n; });
  sumC.foreach ([&] (int n) { *sumC += n; });
  sumG.foreach ([&] (int n) { *sumG += n; });
  sumT.foreach ([&] (int n) { *sumT += n; });
  sumN.foreach ([&] (int n) { *sumN += n; });

  // Dump results
  cout << "|A|=" << *sumA << endl;
  cout << "|C|=" << *sumC << endl;
  cout << "|G|=" << *sumG << endl;
  cout << "|T|=" << *sumT << endl;
  cout << "|N|=" << *sumN << endl;
  
  // Dump some dispatcher values
  cout << "nbCores=" << status.nbCores << "  time=" << status.time << " ms" << endl;
}
//! [snippet1]
