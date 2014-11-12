//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <fstream>
#include <string.h>

/********************************************************************************/
/*                              Minimizers                                      */
/*                                                                              */
/* This snippet shows how to iterate the kmers and have minimizer statistics.   */
/*                                                                              */
/********************************************************************************/

// We use the required packages
using namespace std;

/** Kmer span definition. */
const size_t span = KSIZE_1;

/** Some shortcuts. */
typedef Kmer<span>::ModelDirect                     ModelDirect;
typedef Kmer<span>::ModelCanonical                  ModelCanonical;
typedef Kmer<span>::ModelMinimizer<ModelCanonical>  ModelMinimizer;

typedef ModelMinimizer  Model;

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "you must provide at least 3 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"                   << endl;
        cerr << "   2) minimizer size"              << endl;
        cerr << "   3) FASTA file"                  << endl;
        cerr << "   4) verbose (0/1, 0 by default)" << endl;
        return EXIT_FAILURE;
    }

    // We get the kmer and minimizer sizes.
    size_t kmerSize = atoi(argv[1]);
    size_t mmerSize = atoi(argv[2]);

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        u_int64_t nbKmers       = 0;
        bool display = argc>=5 ? atoi(argv[4]) : false;

        // We declare a Bank instance defined by a list of filenames
        IBank* bank = Bank::singleton().createBank (argv[3]);
        LOCAL (bank);

        // We declare a kmer model and a minimizer model
        Model model (kmerSize, mmerSize);

        // We get a reference on the minimizer model, which will be useful for dumping
       const ModelMinimizer::Model& modelMinimizer = model.getMmersModel();

        Kmer<span>::Type checksum;
        size_t nbChanged = 0;
        size_t nbInvalid = 0;

        // We loop over sequences.
        bank->iterate ([&] (Sequence& seq)
        {
            // We iterate the kmers (and minimizers) of the current sequence.
            model.iterate (seq.getData(), [&] (const Model::Kmer& kmer, size_t idx)
            {
                nbKmers ++;
                if (kmer.hasChanged() == true)   { nbChanged++;  }
                if (kmer.isValid()    == false)  { nbInvalid++;  }
                checksum += kmer.minimizer().value();
            });
        });

        cout << "nbKmers   : " << nbKmers << endl;
        cout << "nbInvalid : " << nbInvalid << endl;
        cout << "nbChanged : " << nbChanged << endl;
        cout << "ratio     : " << (nbChanged > 0 ? (double)nbKmers / (double)nbChanged : 0) << endl;
        cout << "checksum  : " << checksum << endl;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
