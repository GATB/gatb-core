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
/* Cmd-line: kmer5 -in <fasta/q file> -kmer-size <value> ...                    */
/*                  -minimizer-size <size>                                      */
/*                                                                              */
/* Sample: kmer5 -in gatb-core/gatb-core/test/db/reads1.fa -kmer-size 11 ...    */
/*               -minimizer-size 11                                             */
/*                                                                              */
/********************************************************************************/

// We use the required packages
using namespace std;

/** Kmer span definition. */
const size_t span = KMER_SPAN(0);

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
        optimum.setVal(0);
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

typedef ModelMinimizer  Model;

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

        // We get the kmer and minimizer sizes.
        size_t kmerSize = options->getInt(STR_KMER_SIZE);
        size_t mmerSize = options->getInt(STR_MINIMIZER_SIZE);

        // We define a try/catch block in case some method fails (bad filename for instance)
        u_int64_t nbKmers       = 0;
        bool display = options->get(STR_VERBOSE) != 0;

        // We declare a Bank instance defined by a list of filenames
        IBank* bank = Bank::open (options->getStr(STR_URI_INPUT));
        LOCAL (bank);

        // We declare a kmer model and a minimizer model
        Model model (kmerSize, mmerSize);

        // We get a reference on the minimizer model, which will be useful for dumping
       const ModelMinimizer::Model& modelMinimizer = model.getMmersModel();

        Kmer<span>::Type checksum;
        size_t nbChanged = 0;
        size_t nbInvalid = 0;

        // We define an iterator that encapsulates the sequences iterator with progress feedback
        ProgressIterator<Sequence> iter (*bank, "iterate bank");

        // We loop over sequences.
        for (iter.first(); !iter.isDone(); iter.next())
        {
            // Shortcut
            Sequence& seq = iter.item();

//! [snippet1_iterate]
            // We iterate the kmers (and minimizers) of the current sequence.
            model.iterate (seq.getData(), [&] (const Model::Kmer& kmer, size_t idx)
            {
                nbKmers ++;
                if (kmer.hasChanged() == true)   { nbChanged++;  }
                if (kmer.isValid()    == false)  { nbInvalid++;  }
                checksum += kmer.minimizer().value();
            });
//! [snippet1_iterate]
        }

        cout << "nbKmers   : " << nbKmers   << endl;
        cout << "nbInvalid : " << nbInvalid << endl;
        cout << "nbChanged : " << nbChanged << endl;
        cout << "ratio     : " << (nbChanged > 0 ? (double)nbKmers / (double)nbChanged : 0) << endl;
        cout << "checksum  : " << checksum  << endl;
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
