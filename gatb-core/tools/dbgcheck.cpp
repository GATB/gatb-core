#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
struct Info
{
    Info() : checksumNodes(0), checksumSuccessors(0), nbSuccessors(0), abundance(0) {}
    NativeInt64 checksumNodes;
    NativeInt64 checksumSuccessors;
    u_int64_t   nbSuccessors;
    u_int64_t   abundance;
};

/********************************************************************************/
class Functor : public IteratorFunctor
{
public:

    Functor (Graph<NativeInt64>& graph, Info& globalInfo)  : _graph(graph), _globalInfo(globalInfo), _nodes(graph)  {}

    ~Functor ()
    {
        _globalInfo.checksumNodes      += _info.checksumNodes;
        _globalInfo.checksumSuccessors += _info.checksumSuccessors;
        _globalInfo.nbSuccessors       += _info.nbSuccessors;
        _globalInfo.abundance          += _info.abundance;
    }

    void operator ()  (const Node<NativeInt64>& node)
    {
        _info.checksumNodes += node.kmer.value;
        _info.abundance     += node.kmer.abundance;

        /** We retrieve the successors. */
        size_t nbSuccessors = _graph.getSuccessors (node, _nodes);

        _info.nbSuccessors += nbSuccessors;

        /** We iterate all the successors. */
        for (size_t i=0; i<nbSuccessors; i++)   {  _info.checksumSuccessors += _nodes[i].kmer.value;  }
    }

    Graph<NativeInt64>&  _graph;
    NodeSet<NativeInt64> _nodes;
    Info&                _globalInfo;
    Info                 _info;
};

/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser;
    parser.add (new OptionOneParam (STR_URI_INPUT,  "reads input", true));

    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** We load the graph from the provided uri. */
        Graph<NativeInt64> graph = GraphFactory::load <NativeInt64> (options->getStr(STR_URI_INPUT));

        /** We get an iterator over all the nodes of the graph. */
        NodeIterator<NativeInt64> itNodes = graph.nodes();

        /** We encapsulate the nodes iterator with a notification iterator. */
        SubjectIterator<Node<NativeInt64> > itNotif (
            &itNodes,
            itNodes.getNbItems()/100,
            new ProgressTimer(itNodes.getNbItems(), "computing")
        );

        /** We want to compute some statistics about the graph nodes. */
        Info info;

        TIME_START (ti, "compute");

            /** We iterate the nodes through a dispatcher. */
            ParallelCommandDispatcher().iterate (itNotif, Functor(graph, info));

        TIME_STOP (ti, "compute");

        /** We dump the checksum. */
        cout << "nbSuccessors="       << info.nbSuccessors           << "  "
             << "abundance="          << info.abundance              << "  "
             << "checkumNodes="       << info.checksumNodes          << "  "
             << "checksumSuccessors=" << info.checksumSuccessors     << "  "
             << "time="               << ti.getEntryByKey("compute") << "  "
             << endl;
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
    }


}
