// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <queue>
#include <stack>
#include <map>
using namespace std;

#define DEBUG(a)  //a
#define INFO(a)   //a

/********************************************************************************/
class CustomTool : public Tool
{
public:

    // Constructor
    CustomTool () : Tool ("DotGenerator")
    {
        getParser()->push_front (new OptionOneParam (STR_URI_INPUT,  "graph file", true ));
        getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT, "dot file", true ));
    }

    // Actual job done by the tool is here
    void execute ()
    {
        FILE* output = fopen (getInput()->getStr(STR_URI_OUTPUT).c_str(), "w");
        if (output != NULL)
        {
            fprintf (output, "digraph branching {\n");

            // We load the graph
            Graph graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

            map<Node::Value, u_int64_t> mapping;
            u_int64_t count = 0;
            Graph::Iterator<BranchingNode> it = graph.iterator<BranchingNode> ();
            for (it.first(); !it.isDone(); it.next())  { mapping[it.item().kmer] = count++; }

            for (it.first(); !it.isDone(); it.next())
            {
                BranchingNode current = it.item();

                Graph::Vector<BranchingNode> neighbors = graph.successors<BranchingNode> (current);

                for (size_t i=0; i<neighbors.size(); i++)
                {
                    fprintf (output, "%ld -> %ld;\n", mapping[current.kmer], mapping[neighbors[i].kmer]);
                }
            }

            fprintf (output, "}\n");

            fclose (output);
        }

    }
};

/********************************************************************************/
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.
        CustomTool().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
