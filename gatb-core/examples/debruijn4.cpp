//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
struct Info
{
    Info() : checksumNodes(0), checksumSuccessors(0), nbSuccessors(0) {}
    NativeInt64 checksumNodes;
    NativeInt64 checksumSuccessors;
    u_int64_t   nbSuccessors;
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
    }

    void operator ()  (const Node<NativeInt64>& node)
    {
        _info.checksumNodes += node.kmer;

        /** We retrieve the successors. */
        size_t nbSuccessors = _graph.getSuccessors (node, _nodes);

        _info.nbSuccessors += nbSuccessors;

        /** We iterate all the successors. */
        for (size_t i=0; i<nbSuccessors; i++)   {  _info.checksumSuccessors += _nodes[i].kmer;  }
    }

    Graph<NativeInt64>&  _graph;
    NodeSet<NativeInt64> _nodes;
    Info&                _globalInfo;
    Info                 _info;
};

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide:" << endl;
        cerr << "   1) reads file"  << endl;
        cerr << "   2) kmer size"   << endl;
        return EXIT_FAILURE;
    }

    char*  bankUri  = argv[1];
    size_t kmerSize = atoi (argv[2]);
    size_t nks      = argc >= 4 ? atoi (argv[3]) : 3;

    /** We create the graph with  1) a FASTA bank   2) a kmer size */
    Graph<NativeInt64> graph = GraphFactory::createGraph <NativeInt64> (new Bank (bankUri), kmerSize, nks);

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
    cout << "nbNodes="            << info.nbSuccessors           << "  "
         << "checkumNodes="       << info.checksumNodes          << "  "
         << "checksumSuccessors=" << info.checksumSuccessors     << "  "
         << "time="               << ti.getEntryByKey("compute") << "  "
         << endl;

    /** We dump some information about the graph. */
    cout << graph.getInfo() << endl;
}
//! [snippet1]
