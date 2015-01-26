//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <fstream>
#include <queue>
#include <stack>
#include <map>
using namespace std;

#define DEBUG(a)  //a
#define INFO(a)   //a

/********************************************************************************/

const char* STR_NODE_TYPE = "-type";

/********************************************************************************/
class DotGeneratorTool : public Tool
{
public:

    // Constructor
    DotGeneratorTool () : Tool ("DotGenerator")
    {
        _parser->push_front (new OptionOneParam (STR_URI_GRAPH,  "graph file", true ));
        _parser->push_front (new OptionOneParam (STR_URI_OUTPUT, "dot file",  false ));
        _parser->push_front (new OptionOneParam (STR_NODE_TYPE,  "node type (0: all,  1:branching)", false, "1" ));
    }

    template<typename NodeType>
    void process (const char* name)
    {
        string outputFile = getInput()->get(STR_URI_OUTPUT) ?
            getInput()->getStr(STR_URI_OUTPUT) :
            (System::file().getBaseName(getInput()->getStr(STR_URI_GRAPH)) + ".dot");

        ofstream output (outputFile.c_str());

        output << "digraph " << name << "{\n";

        // We load the graph
        Graph graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));

        map<Node, u_int64_t> mapping;
        u_int64_t count = 0;

        Graph::Iterator<NodeType> itMap = graph.iterator<NodeType> ();
        for (itMap.first(); !itMap.isDone(); itMap.next())  { mapping[itMap.item()] = count++; }

        ProgressGraphIterator<NodeType,ProgressTimer> it = graph.iterator<NodeType> ();
        for (it.first(); !it.isDone(); it.next())
        {
            NodeType current = it.item();

            Graph::Vector<NodeType> neighbors = graph.neighbors<NodeType> (current.kmer);

            for (size_t i=0; i<neighbors.size(); i++)
            {
                output << mapping[current.kmer] << " -> " <<  mapping[neighbors[i].kmer] << " ;\n";
            }
        }

        output << "}\n";

        output.close();
    }

    // Actual job done by the tool is here
    void execute ()
    {
        switch (getInput()->getInt(STR_NODE_TYPE))
        {
            case 0: process<Node>          ("all");        break;
            case 1: process<BranchingNode> ("branching");  break;
            default: break;
        }
     }
};

/********************************************************************************/
/*                                                                              */
/*                   Generate dot file from a graph.                            */
/*                                                                              */
/*  This snippet generates a dot file from a graph file. You can then generate  */
/*  a pdf file with "dot -Tpdf graph.dot -o graph.pdf"                          */
/*                                                                              */
/*  NOTE: de Bruijn graphs may be huge and complex, so dot is not the best tool */
/*  to display such graphs. You should use it on small graphs with only a few   */
/*  hundreds of nodes.                                                          */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.
        DotGeneratorTool().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
//! [snippet1]
