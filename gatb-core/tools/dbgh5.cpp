#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/

template <typename T>
class GraphGenerator
{
public:

    static void build (IProperties* options)
    {
        LOCAL (options);

        /** We create the graph. */
        Graph<T> graph = GraphFactory::createGraph <T> (options);

        /** We build the graph. */
        graph.build (new Bank (options->getStr(STR_URI_INPUT)));

        /** We dump some information about the graph. */
        if (options->get(STR_VERBOSE) != 0)  {  cout << graph.getInfo() << endl;  }
    }
};

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.add (new OptionOneParam (STR_URI_INPUT,  "reads input",              true));
    parser.add (new OptionOneParam (STR_URI_OUTPUT, "graph output",             false));
    parser.add (new OptionOneParam (STR_KMER_SIZE,  "kmer size",                false, "27"));
    parser.add (new OptionOneParam (STR_NKS,        "kmer abundance threshold", false, "3" ));
    parser.add (new OptionOneParam (STR_MAX_MEMORY, "max memory",               false, "1000"));
    parser.add (new OptionOneParam (STR_MAX_DISK,   "max disk",                 false, "0"));
    parser.add (new OptionOneParam (STR_NB_CORES,   "nb cores (0 for all)",     false, "0"));
    parser.add (new OptionNoParam  (STR_VERBOSE,    "verbosity",                false));
    parser.add (new OptionNoParam  (STR_HELP,       "help",                     false));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** Shortcut */
        size_t kmerSize = options->getInt (STR_KMER_SIZE);

        /** According to the kmer size, we instantiate one DSKAlgorithm class and delegate the actual job to it. */
        if (kmerSize < 32)        {  GraphGenerator<NativeInt64>::build (options); }
        else if (kmerSize < 64)
        {
    #ifdef INT128_FOUND
            GraphGenerator<NativeInt128>::build (options);
    #else
            GraphGenerator<LargeInt<2> >::build (options);
    #endif
        }
        else if (kmerSize < 96)   {  GraphGenerator<LargeInt<3> >::build (options);  }
        else if (kmerSize < 128)  {  GraphGenerator<LargeInt<4> >::build (options);  }
        else  { throw Exception ("unsupported kmer size %d", kmerSize);  }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
    }
    catch (...)
    {

    }
}
