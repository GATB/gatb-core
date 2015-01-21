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

    void include_frequency (uint32_t* freq_order) {}
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

    void include_frequency (uint32_t* freq_order) {}
};

/** Some shortcuts. */
typedef Kmer<span>::ModelDirect    ModelDirect;
typedef Kmer<span>::ModelMinimizer<ModelDirect,CustomMinimizer> ModelMinimizer;
typedef Kmer<span>::ModelMinimizer<ModelDirect,CustomMaximizer> ModelMaximizer;


/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("KmerTest");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,      "bank input",     true));
    parser.push_back (new OptionOneParam (STR_KMER_SIZE,      "kmer size",      true));
    parser.push_back (new OptionOneParam (STR_MINIMIZER_SIZE, "minimizer size", true));
    parser.push_back (new OptionNoParam  (STR_VERBOSE,        "display kmers",  false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        string bankFilename = options->getStr(STR_URI_INPUT);

        // We get the kmer and minimizer sizes.
        size_t kmerSize = options->getInt(STR_KMER_SIZE);
        size_t mmerSize = options->getInt(STR_MINIMIZER_SIZE);

        // We define a try/catch block in case some method fails (bad filename for instance)
        u_int64_t nbSequences   = 0;
        u_int64_t nbKmers       = 0;
        u_int64_t nbMinimizers  = 0;
        u_int64_t nbMinimizers2 = 0;
        bool display = options->get(STR_VERBOSE) != 0;

        // We declare a bank instance defined by a list of filenames
        IBank* bank = Bank::open (bankFilename);

        // We declare a kmer model and a minimizer model
        ModelMinimizer model (kmerSize, mmerSize);

        // We get a reference on the minimizer model, which will be useful for dumping
        const ModelDirect& modelMinimizer = model.getMmersModel();

        // We create an iterator over this bank.
        ProgressIterator<Sequence> itSeq (*bank);

        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            size_t nbMinimizersPerRead = 0;

            // We iterate the kmers (and minimizers) of the current sequence.
            model.iterate (itSeq->getData(), [&] (const ModelMinimizer::Kmer& kmer, size_t idx)
            {
                // We could display some information about the kmer and its minimizer
                if (display)
                {
                    std::cout << "KMER=" << model.toString(kmer.value()) << "  "
                          << (kmer.hasChanged() ? "NEW" : "OLD") << " "
                          << "MINIMIZER=" << modelMinimizer.toString(kmer.minimizer().value()) << " "
                          << "at position " << kmer.position()
                          << std::endl;
                }

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
        info.add (1, "bank",          "%s",   bankFilename.c_str());
        info.add (1, "kmer_size",     "%ld",  kmerSize);
        info.add (1, "mmer_size",     "%ld",  mmerSize);
        info.add (1, "nb_sequences",  "%ld",  nbSequences);
        info.add (1, "nb_kmers",      "%ld",  nbKmers);
        info.add (1, "nb_minimizers", "%ld",  nbMinimizers);
        info.add (1, "mean",          "%.2f", mean);
        info.add (1, "deviation",     "%.2f", devia);
        cout << info << endl;
    }
    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
