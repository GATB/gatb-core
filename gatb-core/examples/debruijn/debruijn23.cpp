//! [snippet1]

// We include what we need for the test
#include <gatb/gatb_core.hpp>

#include <queue>
#include <stack>
#include <map>
using namespace std;

#define DEBUG(a)  //a
#define INFO(a)   //a

//#undef NDEBUG
#include <cassert>
#define ASSERT(a)

/********************************************************************************/
/********************************************************************************/
#if __SSE3__
    #include <pmmintrin.h>
#else
#if __SSE2__
    #include <emmintrin.h>
#else
    #warning error undefined __SSE3__ or __SSE2__
#endif
#endif

static const size_t PREC = 16;

#if 1
    #define BloomGroupType  BloomGroup
#else
    #define BloomGroupType  BloomGroupCacheCoherent
#endif

#define Counter CounterSSE

/********************************************************************************/

template<size_t PREC>
class CounterBase
{
public:

    CounterBase () : _nbParts(0),  _nbOccurs(0)  {}

    CounterBase (size_t nbParts) : _nbParts(nbParts),  _nbOccurs(0)
    {
        _nbOccurs = new size_t [_nbParts];
    }

    CounterBase (const CounterBase& c) : _nbParts(c._nbParts),  _nbOccurs(0)
    {
        _nbOccurs = new size_t [_nbParts];
    }

    ~CounterBase()  {  if (_nbOccurs) { delete[] _nbOccurs; }  }

protected:
    size_t  _nbParts;
    size_t* _nbOccurs;
};

/********************************************************************************/
/********************************************************************************/
template<size_t PREC>
class CounterSSE : public CounterBase<PREC>
{
public:

    typedef typename BloomGroupType <Kmer<>::Type,PREC>::Result BloomResult;

    CounterSSE ()  {}
    CounterSSE (size_t nbParts) : CounterBase<PREC> (nbParts) {}
    CounterSSE (const CounterBase<PREC>& c) : CounterBase<PREC>(c) {}

    size_t* count (std::vector<BloomResult>& occurs, size_t nbKmers)
    {
        /** We round to the nearest integer modulo 16 (SSE intent...) */
        nbKmers = (nbKmers + 15) & ~15;

        /** We reset the "nb occurs per bloom" table. */
        memset (this->_nbOccurs, 0, sizeof(size_t)*this->_nbParts);

        size_t jmax = this->_nbParts/8;

        /** We want to use a SSE register; in particular we used an union to ease our job,
         * either seeing the register as itself, or as an array of 16 bytes. */
        union { __m128i x; uint8_t b[16]; } v;

        /** The first loop (say X axis) is on the successive kmers of the current sequence. The idea is to
         * vectorize the job for 16 kmers in the same time. */
        for (size_t X=0; X<nbKmers; X+=16)
        {
            /** The second loop (say Y axis) is on the bloom filters. What we got by requesting the Bloom
             * group is an array of bits telling (for one kmer) the presence status in each N Bloom filters
             * for the considered kmer.
             */
            for (size_t Y=0; Y<jmax; Y++)
            {
                /** We fill our SSE register. */
                for (size_t k=0; k<16; k++)
                {
                    u_int8_t* bloomValue = (u_int8_t*) occurs[X+k].array();

                    v.b[k] = bloomValue[Y];
                }

                /** Now we do some computation on this SSE register. */
                for (int l=8; -- l>=0; v.x = _mm_slli_epi16 (v.x, 1))
                {
                    /** We get a 16 bits value (1 for present, 0 otherwise) for the bloom number (8*j+l). */
                    int res = _mm_movemask_epi8 (v.x);

                    /** We count how many bits we got for this bloom filter and increase its kmers occurrence nb. */
                    this->_nbOccurs[8*Y+l] += __builtin_popcount (res);

                    /** Note the last part of the for loop: we shift the 16 integers to 1 bit to the left
                     *  => we will so access to the information of the bloom filter whose index is just below the
                     *  one we just have dealt with. */
                }
            }
        }

        return  this->_nbOccurs;
    }
};

/********************************************************************************/
class Marker
{
public:

    Marker (const Graph& graph) : graph(graph)  { }

    void clear ()  { markMap.clear(); }

    void mark (const Edge& item)
    {
        ASSERT (graph.isBranching(item.from));
        this->markMap[item.from.kmer] |= getMask (item);
    }

    bool isMarked (const Edge& item)
    {
        ASSERT (graph.isBranching(item.from));
        return this->markMap[item.from.kmer]  &  getMask(item);
    }

private:

    const Graph& graph;
    map<Node::Value,u_int8_t>  markMap;

    u_int8_t getMask (const Edge& edge) const  { return (1<<edge.nt);  }
};

/********************************************************************************/
class DFS
{
public:

    DFS (const Graph& graph)  : graph(graph), marker(graph), nbBranching(0) {}

    template<typename Functor>
    void execute (Node node, Functor functor)
    {
        execute_aux (node, functor);
        execute_aux (graph.reverse(node), functor);
    }

    //
    void clear ()  { marker.clear(); }

private:
    const Graph&    graph;
    Marker          marker;
    size_t nbBranching;


    template<typename Functor>
    void execute_aux (Node node, Functor functor)
    {
        stack<Edge> queue;

        // We init the queue with the successors of the provided node
        graph.successors<Edge>(node).iterate ([&] (const Edge& edge)  { queue.push (edge); });

        Graph::Vector<Edge> successors;
        Graph::Vector<Edge> predecessors;
        bool isSimple = true;

        size_t             nbNodes = 0;
        vector<Nucleotide> path (1000);

        while (queue.empty()==false)
        {
            Edge edge        = queue.top();
            Edge currentEdge = edge;
            queue.pop();

            path.clear();

            if (marker.isMarked (edge) == false)
            {
                do  {
                    // We call the functor for the current edge.
                    functor (edge);

                    // We get the successors of the current node.
                    successors =  graph.successors<Edge>(edge.to);

                    // We may get the predecessors of the current node.
                    if (successors.size() == 1) {  predecessors =  graph.predecessors<Edge>(edge.to); }

                    isSimple = (successors.size()==1 && predecessors.size()==1);

                    if (isSimple == true)
                    {
                        path.push_back (edge.nt);
                        edge = successors[0];
                    }

                }  while (isSimple == true);

                marker.mark (currentEdge);
                marker.mark (graph.reverse(edge));

                for (size_t i=0; i<successors.size(); i++)  {  queue.push (successors[i]);  }
            }
        }
    }
};

/********************************************************************************/
class CustomTool : public Tool
{
public:

    // Constructor
    CustomTool () : Tool ("CoverageEstimator")
    {
        getParser()->push_front (new OptionOneParam (STR_URI_INPUT, "graph file", true ));
    }

    // Actual job done by the tool is here
    void execute ()
    {
        // We load the graph
        Graph graph = Graph::load (getInput()->getStr(STR_URI_INPUT));

        // We check that the sorting count got all the kmers
//        if (graph.getInfo().getInt("nks") != 1)  { throw Exception("min abundance must be 1"); }

        // We get an iterator for all nodes of the graph. We use a progress iterator to get some progress feedback
        ProgressGraphIterator<BranchingNode,ProgressTimer>  it (graph.iterator<BranchingNode>(), "graph size ");

        // We create a Depth First Search object
        DFS dfs (graph);

        // We want to count the aggregated size of the unitigs
        u_int64_t unitigsSize = 0;

        // We create a kmer model for iterating kmers of sequences.
        Kmer<>::ModelCanonical model (graph.getKmerSize());

size_t ok=0;
size_t nbBNodes=0;


        u_int64_t nbParts = PREC*64;

        cout << "CREATING BANKS" << endl;
        // We create the N partition bank
        vector<IBank*> subBanks (nbParts);
        for (size_t i=0; i<nbParts; i++)
        {
            stringstream ss;  ss << "part_" << i << ".fa";
            subBanks[i] = BankRegistery::singleton().createBank (ss.str());
        }
        cout << "BANKS CREATED" << endl;


        ////////////////////////////////////////////////////////////
        //
        ////////////////////////////////////////////////////////////

#if 1
        // We loop the branching nodes
        it.iterate ([&] (BranchingNode& node)
        {
            dfs.execute (node, [&] (const Edge& edge)
            {
                unitigsSize++;
            });
        });
#else
        unitigsSize = graph.getInfo().getInt("kmers_nb_solid");
#endif

        cout << "unitigsSize=" << unitigsSize << endl;

        ////////////////////////////////////////////////////////////
        //
        ////////////////////////////////////////////////////////////

        u_int64_t kmersPerPart = unitigsSize / nbParts;

cout << "nbParts=" << nbParts << "  kmersPerPart=" << kmersPerPart << endl;

        /** We define the size of the Bloom filters.*/
        size_t partSize = unitigsSize / nbParts;
        double lg2 = log(2);
        float     NBITS_PER_KMER = log (16*graph.getKmerSize()*(lg2*lg2))/(lg2*lg2);
        size_t    nbHash         = (int)floorf (0.7*NBITS_PER_KMER);
        u_int64_t bloomSize      = (u_int64_t) (partSize * NBITS_PER_KMER);

cout << "NBITS_PER_KMER=" << NBITS_PER_KMER << "  nbHash=" << nbHash << "  bloomSize=" << bloomSize << endl;

        // We clear the dfs
        DFS dfs2 (graph);

        // We loop the branching nodes
        u_int64_t nbNodes = 0;
        size_t    bankIdx = 0;

        u_int64_t maxMemory = (u_int64_t)(4*1000) * (u_int64_t) (1000*1000);
        BloomGroupType<Kmer<>::Type, PREC> bloomGroup (bloomSize, maxMemory, nbHash);

        ProgressGraphIterator<BranchingNode,ProgressTimer>  it2 (graph.iterator<BranchingNode>(), "split graph");
        it2.iterate ([&] (const BranchingNode& node)
        {
            dfs2.execute (node, [&] (const Edge& edge)
            {
                bloomGroup.insert (edge.from.kmer.get<Kmer<>::Type>(), bankIdx);

                nbNodes ++;

                if (nbNodes >= kmersPerPart)
                {
                    bankIdx ++;
                    nbNodes = 0;
                }
            });
        });

        ////////////////////////////////////////////////////////////
        //
        ////////////////////////////////////////////////////////////

        // We retrieve a handle of the bank (we get the uri from the graph properties)
        IBank* bank = BankRegistery::singleton().createBank (graph.getInfo().getStr("input"));
        LOCAL (bank);

        // We create a sequences iterator
        Iterator<Sequence>* iter = this->createIterator<Sequence> (*bank, "iterate bank");
        LOCAL (iter);

        typedef typename BloomGroupType <Kmer<>::Type,PREC>::Result BloomResult;

        std::vector<BloomResult> occurs;

        Counter<PREC> counter (nbParts);

        size_t* nbOccurs = new size_t [nbParts];

        vector<size_t> sizePerPartition (nbParts);

        size_t nbSeqTotal = 0;
        size_t nbSeqOk    = 0;

        u_int64_t nbDataTotal = 0;
        u_int64_t nbDataOk    = 0;

        size_t branchingTheshold = 50;

        // We iterate the bank
        //IDispatcher::Status status = getDispatcher()->iterate (iter, [&] (Sequence& seq)
        iter->iterate ([&] (Sequence& seq)
        {
            nbSeqTotal++;
            nbDataTotal += seq.getDataSize();

            /** We resize the vector of bloom results per kmer of the sequence.
             * Note: we resize it to the nearest integer modulo 16 (SSE intent...) */
            occurs.resize ( (seq.getDataSize() + 15) & ~15);

            /** We reset the vector. */
            std::fill (occurs.begin(), occurs.end(), 0);

            size_t nbKmer = 0;
            size_t nbBranching = 0;

            // We iterate the kmers of the current sequence
            model.iterate (seq.getData(), [&] (const Kmer<>::ModelCanonical::Kmer& kmer, size_t rank)
            {
                if (graph.isBranching (Node::Value(kmer.value())))  {  nbBranching++;  }

                occurs[nbKmer++] = bloomGroup.contains (kmer.value());
            });

            /** We count sequences that have at least one kmer for the required kmer size. */
            if (nbKmer > 0)
            {
#if 1
                /** We count how many time each kmer of the current sequence appears in each Bloom filter. */
                nbOccurs = counter.count (occurs, nbKmer);
#else
                memset (nbOccurs, 0, sizeof(size_t)*nbParts);

                for (size_t i=0; i<nbParts; i++)
                {
                    for (size_t j=0; j<nbKmer; j++)
                    {
                        nbOccurs[i] += (occurs[j][0] & ((u_int64_t)1<<i)) ? 1 : 0;
                    }
                }
#endif
                size_t localSum=0;
                size_t maxValue=0;
                size_t maxIdx  =0;
                for (size_t i=0; i<nbParts; i++)
                {
                    localSum += nbOccurs[i];
                    if (maxValue<nbOccurs[i]) { maxValue=nbOccurs[i]; maxIdx=i; }
                }
                sizePerPartition[maxIdx] +=  seq.getDataSize();

                // We insert the current sequence in the correct partition
                if (nbBranching < branchingTheshold)
                {
                    subBanks[maxIdx]->insert (seq);
                    nbSeqOk++;
                    nbDataOk += seq.getDataSize();
                }

//                printf ("----------------------------  sequence %d  (len=%d  sum=%d nbBranching=%d idx=%d percent=%.3f)\n",
//                    seq.getIndex(), seq.getDataSize(), localSum, nbBranching,
//                    maxIdx, 100.0*(float)maxValue/(float)localSum
//                );
                //for (size_t i=0; i<nbParts; i++)  {  printf ("%2d ", nbOccurs[i]);  if ((i+1)%64==0) { printf ("\n"); }  }
            }
        });

        // We flush the sub banks
        for (size_t i=0; i<nbParts; i++)
        {
            subBanks[i]->flush();
            delete subBanks[i];
        }

        // We save the Bloom group
        bloomGroup.save ("bloom.bin");


printf ("\n------------------ Size per partition ---------------------\n");
for (size_t i=0; i<nbParts; i++)
{
    printf ("%d  %d\n", i, sizePerPartition[i]);
}

        // We gather some statistics.
        getInfo()->add (1, "stats");
        getInfo()->add (2, "nb_splits",    "%ld",  nbParts);
        getInfo()->add (2, "bank_size",    "%ld",  graph.getInfo().getInt("sequences_size"));
        getInfo()->add (2, "unitigs_size", "%ld",  unitigsSize);
        if (unitigsSize > 0)
        {
            getInfo()->add (2, "coverage",     "%.2f", (float) graph.getInfo().getInt("sequences_size") / (float)unitigsSize);
        }
        getInfo()->add (2, "bloom_size", "%ld",  bloomSize);

        getInfo()->add (2, "branchingTheshold", "%ld",  branchingTheshold);

        getInfo()->add (2, "sequences",  "");
        getInfo()->add (3, "nb_total",  "%ld",  nbSeqTotal);
        getInfo()->add (3, "nb_ok",     "%ld",  nbSeqOk);
        getInfo()->add (3, "percent","%.3f",  100.0 * (float)nbSeqOk / (float)nbSeqTotal);

        getInfo()->add (2, "data",  "");
        getInfo()->add (3, "nb_total",  "%ld",  nbDataTotal);
        getInfo()->add (3, "nb_ok",     "%ld",  nbDataOk);
        getInfo()->add (3, "percent","%.3f",  100.0 * (float)nbDataOk / (float)nbDataTotal);
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

    return EXIT_SUCCESS;}
//! [snippet1]
