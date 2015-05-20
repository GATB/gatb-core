//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how to get kmer information from a raw ascii sequence.    */
/*                                                                              */
/********************************************************************************/

int main (int argc, char* argv[])
{
    // We set the maximal span of the kmers. We use here a constant of gatb-core
    // that gives a default span for kmers up to the second kmer size value
    const size_t span = KMER_SPAN(1);

    // We define some nucleotides sequence.
    const char* seq = argc<2 ? "CTACGAATT" : argv[1];

    std::cout << std::endl << "INITIAL STRING IS '" << seq << "'" << std::endl;

    // We set the kmer size as the length of our sequence.
    size_t kmerSize = strlen(seq);

    if (kmerSize >= span)
    {
        std::cerr << "STRING TOO BIG (" << kmerSize << "),  must be less than " << span << std::endl;
        return EXIT_FAILURE;
    }

    // Once we have defined our span, we define some typedefs
    // Note : such definitions are not mandatory but provides better code readability
    typedef Kmer<span>::ModelDirect                    ModelDirect;
    typedef Kmer<span>::ModelCanonical                 ModelCanonical;
    typedef Kmer<span>::ModelMinimizer<ModelCanonical> ModelMinimizer;

    ////////////////////////////////////////////////////////////
    // FIRST EXAMPLE : model direct
    ////////////////////////////////////////////////////////////
    {
//! [snippet1_direct]
        // We declare a kmer model with kmer size big enough to represent our sequence.
        ModelDirect model (kmerSize);

        // We compute the kmer for a given sequence
        ModelDirect::Kmer kmer = model.codeSeed (seq, Data::ASCII);

        std::cout << std::endl;
        std::cout << "-------------------- DIRECT --------------------" << std::endl;
        std::cout << "kmer  value is: " << kmer.value()                 << std::endl;
        std::cout << "kmer string is: " << model.toString(kmer.value()) << std::endl;
//! [snippet1_direct]
    }

    ////////////////////////////////////////////////////////////
    // SECOND EXAMPLE : model canonical
    ////////////////////////////////////////////////////////////
    {
//! [snippet1_canonical]
        // We declare a kmer model with kmer size big enough to represent our sequence.
        ModelCanonical model (kmerSize);

        // We compute the kmer for a given sequence
        ModelCanonical::Kmer kmer = model.codeSeed (seq, Data::ASCII);

        std::cout << std::endl;
        std::cout << "-------------------- CANONICAL --------------------" << std::endl;
        std::cout << "kmer  value is: " << kmer.value()                    << std::endl;
        std::cout << "kmer string is: " << model.toString(kmer.value())    << std::endl;

        /** With this model, we have extra information. */
        std::cout << "forward value  is: " << kmer.forward()                    << std::endl;
        std::cout << "forward string is: " << model.toString(kmer.forward())    << std::endl;
        std::cout << "revcomp value  is: " << kmer.revcomp()                    << std::endl;
        std::cout << "revcomp string is: " << model.toString(kmer.revcomp())    << std::endl;
        std::cout << "used strand is   : " << toString(kmer.strand())           << std::endl;


//! [snippet1_canonical]
    }

    ////////////////////////////////////////////////////////////
    // THIRD EXAMPLE : model minimizer
    ////////////////////////////////////////////////////////////
    {
//! [snippet1_minimizer]
        // We declare a kmer model with kmer size big enough to represent our sequence.
        // Note that we give a second size, which is the size of the minimizers
        ModelMinimizer model (kmerSize, 8);

        // We get a reference on the minimizer model, which will be useful for dumping
        // string value of a minimizer. Recall that 'model' is a model configured with
        // 'kmerSize' but has also to deal with mmers of size kmerSize/2. The 'getMmersModel'
        // method just provides access to this inner model.
        const ModelCanonical& modelMinimizer = model.getMmersModel();

        // We compute the kmer for a given sequence
        ModelMinimizer::Kmer kmer = model.codeSeed (seq, Data::ASCII);

        std::cout << std::endl;
        std::cout << "-------------------- MINIMIZER --------------------" << std::endl;
        std::cout << "kmer  value is: " << kmer.value()                    << std::endl;
        std::cout << "kmer string is: " << model.toString(kmer.value())    << std::endl;

        // With this model, we have extra information.
        std::cout << "forward value  is: " << kmer.forward()                    << std::endl;
        std::cout << "forward string is: " << model.toString(kmer.forward())    << std::endl;
        std::cout << "revcomp value  is: " << kmer.revcomp()                    << std::endl;
        std::cout << "revcomp string is: " << model.toString(kmer.revcomp())    << std::endl;
        std::cout << "used strand is   : " << toString(kmer.strand())           << std::endl;

        // We can also have information about minimizers.
        // Note :  kmer.minimizer() is of type ModelCanonical, ie the type provided as
        // template argument of the ModelMinimizer class.
        std::cout << "minimizer model size       : " << modelMinimizer.getKmerSize() << std::endl;
        std::cout << "minimizer value is         : " << kmer.minimizer().value()     << std::endl;
        std::cout << "minimizer string is        : " << modelMinimizer.toString(kmer.minimizer().value()) << std::endl;
        std::cout << "minimizer position in kmer : " << kmer.position()   << std::endl;
        std::cout << "minimizer changed          : " << kmer.hasChanged() << std::endl;
//! [snippet1_minimizer]
    }
}
//! [snippet1]
