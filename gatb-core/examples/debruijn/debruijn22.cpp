//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

using namespace std;

/********************************************************************************/
/*            Count number of branching nodes in each read of a bank.           */
/*                                                                              */
/* This snippet maps the sequences of a bank on the de Bruijn graph for this    */
/* bank. The idea is to know, for one sequence, how many branching nodes (BN)   */
/* are present in the sequence.                                                 */
/*                                                                              */
/* As a result, we get a distribution <nb branching nodes per kmer, nb occurs>  */
/*                                                                              */
/* This distribution is saved as a HDF5 file. It is possible to extract data    */
/* from it with the HDF5 tools (such as 'h5dump'). For instance, you can dump   */
/* the distribution with gnuplot with the following line:                       */
/*      h5dump -y -d distrib  output.h5  | grep "^\ *[0-9]" | tr -d " "  | paste - - | gnuplot -p -e 'plot  "-" with lines' */
/*                                                                              */
/* Cmd-line: debruijn22 -graph <h5 file> -out <output file>                     */
/*                                                                              */
/* Sample: debruijn22 -graph gatb-core/gatb-core/test/db/celegans_reads.h5 ...  */
/*                    -out /tmp/out.h5                                          */
/* Note:                                                                        */
/*     - '.h5' file contains the HDF5 formatted representation of a de bruijn   */
/*     graph created from a set of reads.                                       */
/*     - a '.h5' file is created using dbgh5 program provided with GATB-Core.   */
/*     Basic use is as follows:                                                 */
/*        dbgh5 -in <fasta/q file> -out <h5 file>                               */
/*     You can also control kmer-size and kmer abundance, see dbgh5 help.       */
/*                                                                              */
/* WARNING ! THIS SNIPPET SHOWS ALSO HOW TO USE LAMBDA EXPRESSIONS, SO YOU NEED */
/* TO USE A COMPILER THAT SUPPORTS THIS FEATURE.                                */
/*                                                                              */
/********************************************************************************/
class BranchingNodeMapping : public Tool
{
public:

    BranchingNodeMapping () : Tool ("BranchingNodeMapping")
    {
        getParser()->push_back (new OptionOneParam (STR_URI_GRAPH,  "graph file",  true));
        getParser()->push_back (new OptionOneParam (STR_URI_OUTPUT, "output file", false, "output"));
    }

    void execute ()
    {
        try
        {
            // We load the graph with the provided graph uri.
            Graph graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));

            // We check that the sorting count got all the kmers
            if (graph.getInfo().getInt("abundance_min") != 1)  { throw Exception("min abundance must be 1"); }

            // We retrieve a handle of the bank (we get the uri from the graph properties)
            IBank* ibank = Bank::open (graph.getInfo().getStr("bank_uri"));
            LOCAL (ibank);

            // We convert the bank to binary form
            BankConverterAlgorithm converter (ibank, graph.getKmerSize(), "bank.bin");
            converter.execute();
            IBank* bank = converter.getResult();
            LOCAL (bank);

            // We create a kmer model for iterating kmers of sequences.
            Kmer<>::ModelCanonical model (graph.getKmerSize());

            size_t totalNbKmers     = 0;
            size_t totalNbBranching = 0;

            // We create a sequences iterator
            Iterator<Sequence>* iter = this->createIterator<Sequence> (*bank, "iterate bank");
            LOCAL (iter);

            // We need a map for building the distribution.
            ThreadObject<map<int,int> > distrib;

            // We iterate the bank
            IDispatcher::Status status = getDispatcher()->iterate (iter, [&] (Sequence& seq)
            {
                // We get a reference on the local distribution for the current thread
                map<int,int>& localDistrib = distrib();

                // We compute the number of kmers in the current sequence
                int nbKmers = seq.getDataSize() - graph.getKmerSize() + 1;

                if (nbKmers > 0)
                {
                    int nbBranching = 0;

                    // We iterate the kmers of the current sequence
                    model.iterate (seq.getData(), [&] (const Kmer<>::ModelCanonical::Kmer& kmer, size_t rank)
                    {
                        // We count the branching nodes.
                        Node node = Node::Value(kmer.value());
                        if (graph.isBranching (node))  {  nbBranching++;  }
                    });

                    // We increase the (local) distribution for this number of branching nodes per sequence
                    localDistrib[nbBranching] ++;

                    // We also increase the number of seen kmers (synchronization required).
                    __sync_fetch_and_add (&totalNbKmers, nbKmers);
                }
            });

            // We merge the (local) distributions filled by each thread into the final distribution
            distrib.foreach ([&] (map<int,int>& localDistrib)
            {
                for (map<int,int>::iterator it = localDistrib.begin(); it != localDistrib.end(); it++)
                {
                    (*distrib)[it->first] += it->second;
                    totalNbBranching += it->first * it->second;
                }
            });

            // We gather some statistics.
            getInfo()->add (1, "stats");

            getInfo()->add (2, "nb_kmers", "");
            getInfo()->add (3, "total",  "%ld", totalNbKmers);
            getInfo()->add (3, "unique", "%ld", graph.getInfo().getInt("kmers_nb_solid"));
            getInfo()->add (3, "ratio",  "%.3f", (float)totalNbKmers / (float)graph.getInfo().getInt("kmers_nb_solid"));

            getInfo()->add (2, "nb_branching", "");
            getInfo()->add (3, "total", "%ld", totalNbBranching);
            getInfo()->add (3, "unique", "%ld", graph.getInfo().getInt("nb_branching"));
            getInfo()->add (3, "ratio",  "%.3f", (float)totalNbBranching / (float)graph.getInfo().getInt("nb_branching"));

            getInfo()->add (2, "percentage",   "");
            getInfo()->add (3, "total",    "%.3f", 100 * (float)totalNbBranching / (float)totalNbKmers);
            getInfo()->add (3, "unique",   "%.3f", 100 * (float)graph.getInfo().getInt("nb_branching") / (float)graph.getInfo().getInt("kmers_nb_solid"));

            getInfo()->add (1, "exec");
            getInfo()->add (2, "time",     "%.3f", (float)status.time / 1000.0);
            getInfo()->add (2, "nb_cores", "%d", status.nbCores);

            // We create the output file
            Storage* storage = StorageFactory(STORAGE_HDF5).create(getInput()->getStr(STR_URI_OUTPUT), true, false);
            LOCAL (storage);

            // We create the collection in the output storage
            typedef Abundance<NativeInt64, u_int32_t> DistribEntry;
            Collection<DistribEntry>& outDistrib = storage->root().getCollection<DistribEntry>("distrib");

            // We dump the distribution in the output file
            for (map<int,int>::iterator it = (*distrib).begin(); it != (*distrib).end(); it++)
            {
                outDistrib.insert (DistribEntry(it->first, it->second));
            }
            outDistrib.flush();
        }
        catch (Exception& e)
        {
            std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
        }
    }
};

/********************************************************************************/
/*            Count number of branching nodes in each read of a bank.           */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
        // We run the branching node mapper
        BranchingNodeMapping().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}
//! [snippet1]
