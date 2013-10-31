//! [snippet1]
// We include what we need for the test

#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
struct Info
{
    Info() : nbSuccessors(0), abundance(0) {}
    Integer   checksumNodes;
    Integer   checksumSuccessors;
    u_int64_t nbSuccessors;
    u_int64_t abundance;
};

/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "you must provide a graph file" << endl;
        return EXIT_FAILURE;
    }

    /** We load the graph from the provided uri. */
    Graph graph = Graph::load (argv[1]);

    /** We get an iterator over all the nodes of the graph. */
    Graph::Iterator<Node> itNodes = graph.iterator<Node>();

    ThreadObject<Info> infos;

    /** We iterate the nodes through a dispatcher. */
    IDispatcher::Status status = Dispatcher().iterate (itNodes, [&] (const Node& node)
    {
        Info& info = infos();

        info.checksumNodes += node.kmer;
        info.abundance     += node.abundance;

        /** We retrieve the successors. */
        Graph::Vector<Node> successors = graph.successors<Node> (node);

        info.nbSuccessors += successors.size();

        /** We iterate all the successors. */
        for (size_t i=0; i<successors.size(); i++)   {  info.checksumSuccessors += successors[i].kmer;  }
    });

    infos.foreach ([&](const Info& i)
    {
        infos->checksumNodes      += i.checksumNodes;
        infos->checksumSuccessors += i.checksumSuccessors;
        infos->nbSuccessors       += i.nbSuccessors;
        infos->abundance          += i.abundance;
    });

    /** We dump the checksum. */
    cout << "nbSuccessors="       << infos->nbSuccessors        << "  "
         << "abundance="          << infos->abundance           << "  "
         << "checkumNodes="       << infos->checksumNodes       << "  "
         << "checksumSuccessors=" << infos->checksumSuccessors  << "  "
         << "time="               << status.time                << "  "
         << endl;

    /** We dump some information about the graph. */
    cout << graph.getInfo() << endl;
}
//! [snippet1]
