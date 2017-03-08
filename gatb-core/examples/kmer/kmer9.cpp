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
/* Cmd-line: kmer9 -in <fasta/q file>                                           */
/*                                                                              */
/* Sample: kmer9 -in gatb-core/gatb-core/test/db/reads1.fa                      */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We create a command line parser.
    IOptionsParser* parser = SortingCountAlgorithm<>::getOptionsParser();
    LOCAL (parser);

    try
    {
        // We parse the user options.
        IProperties* options = parser->parse (argc, argv);

        // We open the input bank
        IBank* bank = Bank::open (options->getStr(STR_URI_INPUT));

        // We create a SortingCountAlgorithm instance.
        SortingCountAlgorithm<> algo (bank, options);

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
