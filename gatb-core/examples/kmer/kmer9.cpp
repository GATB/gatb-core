//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*                              Sorting count                                   */
/*                                                                              */
/* This snippet shows how to count the kmers of a bank by using the sorting     */
/* count algorithm.                                                             */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a command line parser.
    OptionsParser parser ("SortingCount");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,          "bank input",             true));
    parser.push_back (new OptionOneParam (STR_URI_OUTPUT,         "sorting count output",   true));
    parser.push_back (new OptionOneParam (STR_KMER_SIZE,          "kmer size",              false, "31"));
    parser.push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "abundance min",          false, "3"));

    try
    {
        // We parse the user options.
        IProperties* options = parser.parse (argc, argv);

        // We get the options
        size_t kmerSize = options->getInt(STR_KMER_SIZE);
        size_t nks      = options->getInt(STR_KMER_ABUNDANCE_MIN);

        // We open the input bank
        IBank* bank = Bank::open (options->getStr(STR_URI_INPUT));

        // We create an object for storing the couples [kmer,abundance]
        Storage* storage = StorageFactory(STORAGE_HDF5).create(options->getStr(STR_URI_OUTPUT), true, false);   LOCAL (storage);

        // We create a SortingCountAlgorithm instance.
        SortingCountAlgorithm<> algo (storage, bank, kmerSize, make_pair(nks,0));

        // We launch the algorithm
        algo.execute();

        // We display the stats of the sorting count execution
        cout << *algo.getInfo();
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
