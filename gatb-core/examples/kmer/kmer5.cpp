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

/** Implementation of a functor that defines what is a minimizer according to
 *  the list of mmers within a kmer. We mimic here the default behavior
 *  (ie. minimum of the mmers within a kmer) but we could do anything we want to
 *  and we can more generally talk about 'optimizer' rather than 'minimizer'.
 *
 *  The purpose of such a functor is to select a mmer between M mmers in a
 *  specific kmer :
 *
 *      1) the 'init' method is called to initialize the default 'optimum' before looping
 *         over the mmers
 *
 *      2) the operator() method is called for each mmer with the current optimum
 *         value; we can choose here what kind of optimum we want (minimum for instance)
 *         by updating the 'optimum' value
 */
struct CustomMinimizer
{
    template<class Model>  void init (const Model& model, Kmer<span>::Type& optimum) const
    {
        optimum = model.getKmerMax();
    }

    bool operator() (const Kmer<span>::Type& current, const Kmer<span>::Type& optimum) const
    {
        return current < optimum;
    }
};

/** We define here a 'maximizer' in the mmers of a specific kmer. */
struct CustomMaximizer
{
    template<class Model>  void init (const Model& model, Kmer<span>::Type& optimum) const
    {
        optimum = Kmer<span>::Type(0);
    }

    bool operator() (const Kmer<span>::Type& current, const Kmer<span>::Type& optimum) const
    {
        return !(current < optimum);
    }
};

/** Some shortcuts. */
typedef Kmer<span>::ModelDirect    ModelDirect;
typedef Kmer<span>::ModelMinimizer<ModelDirect,CustomMinimizer> ModelMinimizer;
typedef Kmer<span>::ModelMinimizer<ModelDirect,CustomMaximizer> ModelMaximizer;


/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "you must provide at least 2 arguments. Arguments are:" << endl;
        cerr << "   1) kmer size"           << endl;
        cerr << "   2) minimizer size"      << endl;
        cerr << "   2) list of FASTA files" << endl;
        return EXIT_FAILURE;
    }

    // We get the kmer and minimizer sizes.
    size_t kmerSize = atoi(argv[1]);
    size_t mmerSize = atoi(argv[2]);

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        u_int64_t nbSequences   = 0;
        u_int64_t nbKmers       = 0;
        u_int64_t nbMinimizers  = 0;
        u_int64_t nbMinimizers2 = 0;

        // We declare a Bank instance defined by a list of filenames
        BankFasta b (argc-3, argv+3);

        // We declare a kmer model and a minimizer model
        ModelMinimizer model (kmerSize, mmerSize);

        // We create an iterator over this bank.
        ProgressIterator<Sequence> itSeq (b);

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            size_t nbMinimizersPerRead = 0;

            // We iterate the kmers (and minimizers) of the current sequence.
            model.iterate (itSeq->getData(), [&] (const ModelMinimizer::Kmer& kmer, size_t idx)
            {
                // We may have to update the number of different minimizer in the current sequence
                if (kmer.hasChanged()==true)  {  nbMinimizersPerRead++;  }

                nbKmers ++;
            });

            // We update global statistics
            nbSequences   ++;
            nbMinimizers  += nbMinimizersPerRead;
            nbMinimizers2 += nbMinimizersPerRead*nbMinimizersPerRead;
        }

        double mean  = nbSequences>0 ? (double)nbMinimizers  / (double) nbSequences : 0;
        double mean2 = nbSequences>0 ? (double)nbMinimizers2 / (double) nbSequences : 0;
        double devia = sqrt (mean2 - mean*mean);

        // We dump results
        Properties info;
        info.add (0, "info");
        info.add (1, "bank",          "%s",   argv[3]);
        info.add (1, "kmer_size",     "%ld",  kmerSize);
        info.add (1, "mmer_size",     "%ld",  mmerSize);
        info.add (1, "nb_sequences",  "%ld",  nbSequences);
        info.add (1, "nb_kmers",      "%ld",  nbKmers);
        info.add (1, "nb_minimizers", "%ld",  nbMinimizers);
        info.add (1, "mean",          "%.2f", mean);
        info.add (1, "deviation",     "%.2f", devia);
        cout << info << endl;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
