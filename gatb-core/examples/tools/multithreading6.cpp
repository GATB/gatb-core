//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*                         Multithreaded iteration of a bank.                   */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "you must provide at least the FASTA file path." << endl;
        return EXIT_FAILURE;
    }

    // We get a handle on a bank
    BankFasta bank (argv[1]);

    // We get the number of cores to be used.
    size_t nbCores = (argc >=3 ? atoi(argv[2]) : 0);

    // We create a dispatcher (use all cores by default).
    Dispatcher dispatcher (nbCores);

    // We will count nucleotides occurrences.
    ThreadObject<int> sumA, sumC, sumG, sumT, sumN;

    // We iterate the bank. Note how we provide a bank iterator to the dispatcher
    dispatcher.iterate (bank.iterator(), [&] (const Sequence& seq)
    {
        // We use shortcuts references for the different local sums. It avoids to retrieve
        // them each time a nucleotide of the sequence is handled (see for loop below)
        // and may give much better performance.
        int& localA = sumA();
        int& localC = sumC();
        int& localG = sumG();
        int& localT = sumT();
        int& localN = sumN();

        // We loop the nucleotides of the current sequence.
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

    sumA.foreach ([&] (int n) { *sumA += n; });
    sumC.foreach ([&] (int n) { *sumC += n; });
    sumG.foreach ([&] (int n) { *sumG += n; });
    sumT.foreach ([&] (int n) { *sumT += n; });
    sumN.foreach ([&] (int n) { *sumN += n; });

    cout << "|A|=" << *sumA << endl;
    cout << "|C|=" << *sumC << endl;
    cout << "|G|=" << *sumG << endl;
    cout << "|T|=" << *sumT << endl;
    cout << "|N|=" << *sumN << endl;
}
//! [snippet1]
