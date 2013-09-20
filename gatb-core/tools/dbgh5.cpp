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

        IBank* bank = new Bank (options->getStr(STR_URI_INPUT));

        /** We create the graph. */
        Graph<T> graph = GraphFactory::createGraph <T> (bank, options);

        /** We build the graph. */
        graph.build ();

        /** We dump some information about the graph. */
        cout << graph.getInfo() << endl;
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
