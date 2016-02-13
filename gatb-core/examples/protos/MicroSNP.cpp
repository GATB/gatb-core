
// We include what we need for the test
#include <gatb/gatb_core.hpp>

/********************************************************************************/
class MicroSNP : public Tool
{
public:

    MicroSNP () : Tool ("MicroSNP")
    {
        getParser()->push_back (Graph::getOptionsParser());
    }

    void execute ()
    {
        FILE * snps = fopen ("snps","w");

        // We create the de Bruijn graph
        Graph graph = Graph::create (getParser()->getProperties());

        // We get an iterator for all nodes of the graph. We use a progress iterator to get some progress feedback
        ProgressGraphIterator<BranchingNode,ProgressTimer>  it (graph.iteratorBranching(), "MiniDiscoSNP: finding bubbles          ");

        int nbsnps = 1 ;
        int ksize  = graph.getKmerSize();

        for (it.first(); !it.isDone(); it.next())
        {
            Node& current = it.item();

            int outdegree = graph.outdegree(current);
            int indegree  = graph.indegree (current);

            if ((indegree ==2 && outdegree==1) )
            {
                current = graph.reverse(current);
            }

            if ((indegree ==1 && outdegree==2)  ||  (indegree ==2 && outdegree==1) )
            {
                //get neighbor branching edges
                Graph::Vector<BranchingEdge> branchingNeighbors = graph.successorsBranchingEdge (current);

                //clean bubble
                if(branchingNeighbors.size()==2  && branchingNeighbors[0].distance == ksize && branchingNeighbors[1].distance == ksize)
                {
                    for (size_t i=0; i<branchingNeighbors.size(); i++)
                    {
                        BranchingEdge edge = branchingNeighbors[i];

                        // We avoid duplicates
                        if( edge.to.kmer < edge.from.kmer )  { continue; }

                        flockfile(snps);
                        fprintf(snps,">SNP %i path %zu\n",nbsnps,i );
                        fprintf(snps,"%s%c%s\n",graph.toString (edge.from).c_str(), ascii(edge.nt) , graph.toString (edge.to).c_str() );
                        funlockfile(snps);

                        __sync_fetch_and_add ( &nbsnps, i);
                    }
                }
            }
        }

        getInfo()->add (1, "stats");
        getInfo()->add (2, "nb", "%d", nbsnps);

        fclose(snps);
    }
};

/********************************************************************************/

// Once our tool class is defined, we can run it in the main function of the program.
int main (int argc, char* argv[])
{
    // We use a try/catch block since GATB functions may throw exceptions
    try
    {
        MicroSNP().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}
