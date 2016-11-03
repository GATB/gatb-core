// 
// GATB-Core Online Tutorial 
/********************************************************************************/
/*                    Making a basic SNP detection tool                         */
/*                                                                              */
/* This snippet shows how to locate isolated SNPs from some reads. Actually,    */
/* it illustrates how to implement the basics of DiscoSNP/KisSnp2 algorithm     */
/* as presented in the DiscoSNP publication:                                    */
/*                                                                              */
/*   http://nar.oxfordjournals.org/content/43/2/e11.full#sec-2                  */
/*                                                                              */
/********************************************************************************/

// We include GATB-Core
#include <gatb/gatb_core.hpp>

/********************************************************************************/
// It is worth noting that we use the GATB-Core Tool's API
class MicroSNP : public Tool
{
  public:

  MicroSNP () : Tool ("MicroSNP")
  {
  // If you use the cmd-line, then un-comment the following line
  //getParser()->push_back (Graph::getOptionsParser());
  }

  void execute ()
  {
    FILE * snps = stdout;//fopen ("snps","w");

    // We create the de Bruijn graph.
    // case 1: for the purpose of the online GATB-Tutorial, input file is 
    // provided, so we "force the command-line.
    Graph graph = Graph::create (
      "-in read.fasta -kmer-size 7 -abundance-min 1 -verbose 0");

    // case 2: use the following line when passing arguments from the cmd-line
    // (and also uncomment code line in the constructor)
    //Graph graph = Graph::create (getParser()->getProperties());

    // We get an iterator for all branching nodes of the graph. 
    // We use a progress iterator to get some progress feedback 
    // (will be displayed on stdout)
    ProgressGraphIterator<BranchingNode,ProgressTimer>  it 
      (graph.iteratorBranching(), "MiniDiscoSNP: finding bubbles");

    // used to store #SNPs found
    int nbsnps = 0 ;
  
    // get the kmer size
    int ksize  = graph.getKmerSize();

    // now, we iterate over all branching nodes
    for (it.first(); !it.isDone(); it.next())
    {
      // get the current node into an appropriate variable
      Node& current = it.item();

      // get branching degrees of the node
      int outdegree = graph.outdegree(current);
      int indegree  = graph.indegree (current);

      ///reverse if we have an ending node
      if ( (indegree ==2 && outdegree==1) )
      {
        current = graph.reverse(current);
      }

      // Do we have a possible bubble? Two cases as follows:
      if ((indegree ==1 && outdegree==2)  ||  (indegree ==2 && outdegree==1) )
      {
        // Get neighbor branching edges
        Graph::Vector<BranchingEdge> branchingNeighbors = 
          graph.successorsBranchingEdge (current);

        // Do we have a clean bubble?
        // True if: exactly 2 branches AND size of each branch == kmer size
        if( branchingNeighbors.size()==2 && 
            branchingNeighbors[0].distance == ksize && 
            branchingNeighbors[1].distance == ksize)
        {
          for (size_t i=0; i<branchingNeighbors.size(); i++)
          {
            BranchingEdge edge = branchingNeighbors[i];

            // We avoid duplicates
            if( edge.to.kmer < edge.from.kmer ) { 
              continue;
            }
            // We dump the SNP found
            flockfile(snps);
            fprintf(snps, "\n>SNP %i path %zu\n", nbsnps, i );
            fprintf(snps, " %s%c%s\n", 
              graph.toString(edge.from).c_str(), 
              ascii(edge.nt), 
              graph.toString(edge.to).c_str() );
            funlockfile(snps);
            // We safely increment 'nbsnps'. Remember that GATB Iterators use 
            // multi-threading!
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

// test on this file
//>read_1
//TACACGTCGGCACATCG
//>read_2
//TACACGTCTGCACATCG
//
//with:
//MicroSNP -in read.fasta -kmer-size 7 -abundance-min 1

