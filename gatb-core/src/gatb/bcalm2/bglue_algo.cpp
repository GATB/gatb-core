/* remaining issue:
- no more than 2^(32-1) sequences to glue together (should be ok for spruce)
*/
#include "bglue_algo.hpp"

#include <unordered_map>
#include "unionFind.hpp"
#include <BooPHF/BooPHF.h>
#include "ThreadPool.h"

#include "logging.hpp"
/*#include "buffer_allocator.tcc"
#include "buffer_manager.tcc"*/
#include <sstream>
#include <iomanip>

/*#include "ctpl_stl.h" // alternative to threadpool // https://github.com/vit-vit/CTPL/blob/master/ctpl_stl.h // didn't commit because didnt use
#include "buffer_allocator.h" // memory pool from https://github.com/vincetse/allocator, didn't commit the files because didnt use
#include "buffer_manager.h" // same, memory pool from https://github.com/vincetse/allocator/blob/master/include/lazy/memory/
*/

#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/storage/impl/StorageTools.hpp>

#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/Banks.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/kmer/impl/PartiInfo.hpp>   // for repartitor 
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/collections/impl/BooPHF.hpp>

#include <queue> // for priority_queue


//heh at this point I could have maybe just included gatb_core.hpp but well, no circular dependencies, this file is part of gatb-core now.

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;



using namespace std;


namespace gatb { namespace core { namespace debruijn { namespace impl  {

    template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 1)
{
        std::ostringstream out;
            out << std::fixed << std::setprecision(n) << a_value;
                return out.str();
}

  // a hash wrapper for hashing kmers in Model form
    template <typename ModelType>
    class Hasher_T
    {
       public:
        ModelType model;

        Hasher_T(ModelType &model) : model(model) {};

        // fun fact: I tried with a mask = (1<<25)-1, 
        // and with chr14, it produced one big partition. So i guess that hash image needs to be large
        uint64_t operator ()  (const typename ModelType::Kmer& key, uint64_t seed = 0) const  {
                return model.getHash(key.value()) ;
                }
    };


 template <typename T>
void free_memory_vector(std::vector<T> &vec)
{
    vec.clear();
    vector<T>().swap(vec); // it's a trick to properly free the memory, as clear() doesn't cut it (http://stackoverflow.com/questions/3477715/c-vectorclear)
}


static 
char rc /*cheap desambiguation compared to GraphUnitigs because TemplateSpecialization8 complains */(char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	else if (s == 'a') return 't';
	else if (s == 'c') return 'g';
	else if (s == 'g') return 'c';
	else if (s == 't') return 'a';
	return 'X';
}


static string rc(const string &s) {
	string rcs = "";
	for (signed int i = s.length() - 1; i >= 0; i--) {rcs += rc(((char)s[i]));}
	return rcs;
}


// manipulation of abundance vectors encoded as strings. recent addition (post-publication)

static float get_mean_abundance(const string& list)
{
    float mean_abundance=0;
    int n=0;
    std::stringstream stream(list);
    while(1) {
        int a;
        stream >> a;
        if(!stream)
            break;
        mean_abundance +=a;
        n++;
    }
    return mean_abundance / (float)n;
}

static uint64_t get_sum_abundance(const string& list)
{
    uint64_t sum_abundances=0;
    std::stringstream stream(list);
    while(1) {
        int a;
        stream >> a;
        if(!stream)
            break;
        sum_abundances +=a;
    }
    if (sum_abundances > 2000000000LL) std::cout << "warning, large abundance reached, may have printing problems" << std::endl; // maybe will disrupt optimizing the code of that function, but it's not that critical
    return sum_abundances;
}


static string reverse_abundances(const string& list)
{
    string rev="";
    std::stringstream stream(list);
    while(1) {
        int a;
        stream >> a;
        if(!stream)
            break;
        rev = to_string(a) + " " + rev;
    }
    return rev;
}

static string skip_first_abundance(const string& list)
{
    string res="";
    std::stringstream stream(list);
    bool skip = true;
    while(1) {
        int a;
        stream >> a;
        if(!stream)
            break;
        if (skip)
        {
            skip = false;
            continue;
        }
        res += to_string(a) + " ";
    }
    return res;
}


template<int SPAN>
struct markedSeq
{
    // there used to be "string seq; string abundance" but i noticed that i did not need that info for determining the chain of glues. not much space saved though (like 10-20%). I suppose the biggest memory-hog is the ks/ke unordered_map
    uint64_t index;
    bool rc;
    bool lmark, rmark;
    typedef typename Kmer<SPAN>::Type Type;
    Type ks, ke; // [start,end] kmers of seq, in canonical form (redundant information with seq, but helpful)

    markedSeq(uint64_t index, bool lmark, bool rmark, const Type &ks, const Type &ke) : index(index), rc(false), lmark(lmark), rmark(rmark), ks(ks), ke(ke) {};

    void revcomp()
    {
        rc = !rc;
        std::swap(lmark, rmark);
        std::swap(ks, ke);
    }
};


// hack to refer to a sequences in msInPart as reverse complemented

static uint32_t is_rev_index(uint32_t index)
{
    return ((index >> 31 & 1) == 1);
}


static uint32_t rev_index(uint32_t index)
{
    if (is_rev_index(index))
    { std::cout << "Error: glue sequence index too large " << index << std::endl; exit(1);}
    return index | (1<<31);
}

static uint32_t no_rev_index(uint32_t index)
{
    return index & ((1LL<<31) - 1LL);
}

//typedef lazy::memory::buffer_allocator<markedSeq> custom_allocator_t;
//typedef std::allocator<markedSeq> custom_allocator_t;

/* input: markedSequences, list of sequences in a partition
 * output: res, a list of lists of sequences that will be glued together
 */
template<int SPAN>
static void determine_order_sequences(vector<vector<uint32_t>> &res, const vector<markedSeq<SPAN>> &markedSequences, int kmerSize, bool debug=false)
{
    typedef typename Kmer<SPAN>::Type Type;
    unordered_map<Type, set<uint32_t> > kmerIndex;
    set<uint32_t> usedSeq;
    unsigned int nb_chained = 0;

    // index kmers to their seq
    // kmerIndex associates a kmer extremity to its index in markedSequences
    for (uint32_t i = 0; i < markedSequences.size(); i++)
    {
        kmerIndex[markedSequences[i].ks].insert(i);
        kmerIndex[markedSequences[i].ke].insert(i);
    }

    auto glue_from_extremity = [&](markedSeq<SPAN> current, uint32_t chain_index, uint32_t markedSequence_index)
    {
        vector<uint32_t> chain;
        chain.push_back(chain_index);

        bool rmark = current.rmark;
        usedSeq.insert(markedSequence_index);

        while (rmark)
        {
            if (debug)
                std::cout << "current ke " << current.ke << " index " << no_rev_index(chain_index) << " markings: " << current.lmark << current.rmark <<std::endl;

            // this sequence has a rmark, so necessarily there is another sequence to glue it with. find it here.
            set<uint32_t> candidateSuccessors = kmerIndex[current.ke];
           
            assert(candidateSuccessors.find(markedSequence_index) != candidateSuccessors.end()); // remove the current seq from our indexing data structure 
            candidateSuccessors.erase(markedSequence_index);

            assert(candidateSuccessors.size() == 1); // normally there is exactly one sequence to glue with

            uint32_t successor_index = *candidateSuccessors.begin(); // pop()
            assert(successor_index != markedSequence_index);
            markedSeq<SPAN> successor = markedSequences[successor_index];

            chain_index = markedSequences[successor_index].index;

            if (successor.ks != current.ke || (!successor.lmark))
            {
                successor.revcomp();
                chain_index = rev_index(chain_index);
            }

            // some checks
            {
                if (debug)
                    std::cout << "successor " << successor_index /*<<" successor ks ke "  << successor.ks << " "<< successor.ke*/ /* need to convert Type to string to print, didn't bother writing that code yet */ << " markings: " << successor.lmark << successor.rmark << std::endl;
                assert(successor.lmark);
                assert(successor.ks == current.ke);
                // edge case where the seq to be glued starts and ends with itself. 
                // it should be a kmer (is tested below with an assert())
                if (successor.ks == successor.ke)
                {
                    if (successor.lmark == false)
                        assert(successor.rmark == true);
                    else
                        assert(successor.rmark == false);
                    // it's the only possible cases I can think of
                    // there is actually nothing to be done now, it's an extremity, so it will end.
                    // on a side note, it's pointless to save this kmer in bcalm.
                }
            }

            current = successor;
            markedSequence_index = successor_index;
            chain.push_back(chain_index);
            rmark = current.rmark;
            assert((usedSeq.find(markedSequence_index) == usedSeq.end()));
            usedSeq.insert(markedSequence_index);
        }

        res.push_back(chain);
        nb_chained += chain.size();

    };

    // iterated markedSequences, and picks extremities of a chain
    for (unsigned int i = 0; i < markedSequences.size(); i++)
    {
        markedSeq<SPAN> current = markedSequences[i];
        if (usedSeq.find(i) != usedSeq.end())
        {
            if (debug)
                std::cout << "sequence has already been glued" << std::endl;
            continue; 
        }

        if (current.lmark && current.rmark)
        {
            if (debug)
                std::cout << "not the extremity of a chain" << std::endl;
            continue;  
        }
    
        /* normalize sequence so that lmark is false */
        uint32_t chain_index = markedSequences[i].index;
        if (current.lmark)
        {
            current.revcomp(); 
            chain_index = rev_index(chain_index);
        }

        assert(current.lmark == false);    

        glue_from_extremity(current, chain_index, i);

    }

    /* // attempt to fix circular contigs, but I was in a hurry, so not finished
    while (nb_chained < sequences.size())
    {
        // there might be a special case: a circular unitig, to be glued with multiple sequences all containing doubled kmers at extremities
        // my fix plan: we pick an extremity at random, and chop the last nucleotide and mark it to not be glued. also find the corresponding kmer in other extremity, and mark it as not to be glued

        vector<int> remaining_indices;
        for (uint32_t i = 0; i < markedSequences.size(); i++)
        {
            if (usedSeq.find(i) != usedSeq.end())
                remaining_indices.push_back(i);
        }
        uint32_t chain_index = remaining_indices[0];
        string kmer = markedSequences[chain_index].substr(0,kmerSize);
       // markedSequences[chain_index] = // TODO continue
    }
    */
    if (nb_chained < markedSequences.size())
    {
        std::cout << " Note: " << markedSequences.size() - nb_chained << " sequence chunks not returned in output unitigs (likely small circular contigs)" << std::endl;
    }
    // assert(sequences.size() == nb_chained); // make sure we've scheduled to glue all sequences in this partition
}

/* straightforward glueing of a chain
 * sequences should be ordered and in the right orientation
 * so, it' just a matter of chopping of the first kmer
 */
static void glue_sequences(vector<uint32_t> &chain, std::vector<std::string> &sequences, std::vector<std::string> &abundances, int kmerSize, string &res_seq, string &res_abundances)
{
    bool debug=false;

    string previous_kmer = "";
    unsigned int k = kmerSize;
    
    if (debug) std::cout << "glueing new chain: ";
    for (auto it = chain.begin(); it != chain.end(); it++)
    {
        uint32_t idx = *it;

        string seq = sequences[no_rev_index(idx)];
        string abs = abundances[no_rev_index(idx)];

        if (is_rev_index(idx))
        {
            seq = rc(seq);
            abs = reverse_abundances(abs);
        }
        
        if (previous_kmer.size() == 0) // it's the first element in a chain
        {
            res_seq += seq;
            res_abundances += abs;
        }
        else
        {
            assert(seq.substr(0, k).compare(previous_kmer) == 0);
            res_seq += seq.substr(k);
            res_abundances += skip_first_abundance(abs);
        }
    
        if (debug) std::cout << seq << " ";

        previous_kmer = seq.substr(seq.size() - k);
        assert(previous_kmer.size() == k);
    }
     if (debug) std::cout << std::endl;
}


static void output(const string &seq, gatb::core::debruijn::impl::BufferedFasta &out, const string comment = "")
{
    out.insert(seq, comment);
    // BufferedFasta takes care of the flush
}



 // used to get top N elements of a vector
template <typename T>
struct Comp{
    Comp( const vector<T>& v ) : _v(v) {}
    bool operator ()(T a, T b) { return _v[a] > _v[b]; }
    const vector<T>& _v;
};


// taken from GATB's MPHF.hpp and BooPHF.hpp (except that we don't need the iteration stuff from that file)
template<typename Key>
class hasher_t
{   
    typedef jenkins64_hasher BaseHasher; /* from BooPHF.hpp, which itself is from emphf:base_hasher */
    BaseHasher emphf_hasher;
    AdaptatorDefault<Key> adaptor;

    public:
    hasher_t(){
        std::mt19937_64 rng(37); // deterministic seed
        emphf_hasher = BaseHasher::generate(rng);
    }

    uint64_t operator ()  (const Key& key, uint64_t seed = 0) const  {
        if (seed != 0x33333333CCCCCCCCULL)
            return std::get<0>(emphf_hasher(adaptor(key)));
        return std::get<2>(emphf_hasher(adaptor(key)));
        // this is a big hack, because I'm lazy. 
        // I wanted to return two different hashes depending on how boophf calls it
        // since I contrl BooPHF code's, I know it calls this function with 0x33333333CCCCCCCCULL as the second seed.
    }
};
    
typedef uint64_t partition_t;

/* computes and uniquifies the hashes of marked kmers at extremities of all to-be-glued sequences */
template <int SPAN>
void prepare_uf(std::string prefix, IBank *in, const int nb_threads, int& kmerSize, int pass, int nb_passes, uint64_t &nb_elts, uint64_t estimated_nb_glue_sequences)
{
  
    std::atomic<unsigned long> nb_marked_extremities, nb_unmarked_extremities; 
    nb_marked_extremities = 0; nb_unmarked_extremities = 0;

    std::vector<std::vector<partition_t >> uf_hashes_vectors(nb_threads);
    
    // relatively accurate number of sequences to be inserted
    for (int i = 0; i < nb_threads; i++)
        uf_hashes_vectors[i].reserve(estimated_nb_glue_sequences/(nb_passes*nb_threads));

    /* class (formerly a simple lambda function) to process a kmer and decide which bucket(s) it should go to */
    /* needed to make it a class because i want it to remember its thread index */
    class UniquifyKeys 
    {
        typedef typename Kmer<SPAN>::ModelCanonical ModelCanon;
        
        int k;
        int pass, nb_passes, nb_threads;
        ModelCanon modelCanon;
        Hasher_T<ModelCanon> hasher;
        std::atomic<unsigned long> &nb_marked_extremities, &nb_unmarked_extremities; 
        std::vector<std::vector<partition_t >> &uf_hashes_vectors;
        int _currentThreadIndex;

        public: 
        UniquifyKeys(int k, int pass, int nb_passes, int nb_threads,
                     std::atomic<unsigned long> &nb_marked_extremities, std::atomic<unsigned long> & nb_unmarked_extremities,
                    std::vector<std::vector<partition_t >> &uf_hashes_vectors
                     ) : k(k), pass(pass), nb_passes(nb_passes), nb_threads(nb_threads), modelCanon(k), hasher(modelCanon),
                        nb_marked_extremities(nb_marked_extremities), nb_unmarked_extremities(nb_unmarked_extremities),
                         uf_hashes_vectors(uf_hashes_vectors), _currentThreadIndex(-1)
        {}
 
        void operator()     (const Sequence& sequence) {
            const string seq = sequence.toString(); 
            const string comment = sequence.getComment();

            const bool lmark = comment[0] == '1';
            const bool rmark = comment[1] == '1';
            int thread = getThreadIndex();

            if (lmark)
            {
                const string kmerBegin = seq.substr(0, k );
                // UF of canonical kmers in ModelCanon form, then hashed
                const typename ModelCanon::Kmer kmmerBegin = modelCanon.codeSeed(kmerBegin.c_str(), Data::ASCII);
                const uint64_t h1 = hasher(kmmerBegin);
                if (h1 % (uint64_t)nb_passes == (uint64_t)pass)
                {
                    uf_hashes_vectors[thread].push_back(h1);
                    nb_marked_extremities++;
                }
            }
            else 
                nb_unmarked_extremities++;

            if (rmark)
            {
                const string kmerEnd = seq.substr(seq.size() - k , k );
                const typename ModelCanon::Kmer kmmerEnd = modelCanon.codeSeed(kmerEnd.c_str(), Data::ASCII);
                const uint64_t h2 = hasher(kmmerEnd);

                if (h2 % (uint64_t)nb_passes == (uint64_t)pass)
                {
                    uf_hashes_vectors[thread].push_back(h2);
                    nb_marked_extremities++;
                }
            }
            else
                nb_unmarked_extremities++;
        }

        /* neat trick taken from erwan's later work in gatb to find the thread id of a dispatched function */
        int getThreadIndex()
        {
            if (_currentThreadIndex < 0)
            {
                std::pair<IThread*,size_t> info;
                if (ThreadGroup::findThreadInfo (System::thread().getThreadSelf(), info) == true)
                {
                    _currentThreadIndex = info.second;
                }
                else
                {
                    throw Exception("Unable to find thread index during InsertIntoQueues");
                }
            }
            return _currentThreadIndex;
        }
    };


    Dispatcher dispatcher (nb_threads);

    UniquifyKeys uniquifyKeys(kmerSize, pass, nb_passes, nb_threads,
                              nb_marked_extremities, nb_unmarked_extremities,
                              uf_hashes_vectors);
    dispatcher.iterate (in->iterator(), uniquifyKeys);
    logging( std::to_string(nb_marked_extremities.load()) + " marked kmers, " + std::to_string(nb_unmarked_extremities.load()) + " unmarked kmers");


    //
    // single-threaded version
    /* auto it = in->iterator();    
     for (it->first (); !it->isDone(); it->next())
        prepareUF(it->item());*/


    logging("created vector of hashes, size approx " + std::to_string( sizeof(partition_t)*nb_marked_extremities.load()/1024/1024) + " MB)");
    ThreadPool uf_sort_pool(nb_threads); // ThreadPool
  //  ctpl::thread_pool uf_merge_pool(nb_threads);

    // sort and uniquify UF vectors (from uf_hashes_vector). the uniquify is actually optional but doesn't cost much
    for (int i = 0; i < nb_threads; i++)
    {
        auto sortuniq = [&uf_hashes_vectors, i] (int thread_id)
        {
            std::vector<partition_t> &vec = uf_hashes_vectors[i];
            sort( vec.begin(), vec.end() );
            vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
        };
        uf_sort_pool.enqueue(sortuniq);  // ThreadPool
        //uf_sort_pool.push(sortuniq);  // ctpl
        //sortuniq(0); // single-threaded
    }

    uf_sort_pool.join(); // ThreadPool
    
    
    // a single-threaded merge and write to file before they're loaded again in bglue
    BagFile<uint64_t> * bagf = new BagFile<uint64_t>( prefix+".glue.hashes."+ to_string(pass)); LOCAL(bagf); 
	Bag<uint64_t> * currentbag =  new BagCache<uint64_t> (  bagf, 10000 ); LOCAL(currentbag);// really? we have to through these hoops to do a simple binary file in gatb? gotta change this.
    uint64_t nb_elts_pass = 0;

    priority_queue<std::tuple<uint64_t,int>, std::vector<std::tuple<uint64_t,int>>, std::greater<std::tuple<uint64_t,int>> > pq; // http://stackoverflow.com/questions/2439283/how-can-i-create-min-stl-priority-queue
    vector<uint64_t> hash_vector_idx(nb_threads);
    vector<uint64_t> hash_vector_size(nb_threads);

    // prime the pq
    for (int i = 0; i < nb_threads; i++)
    {
        hash_vector_idx[i] = 0;
        hash_vector_size[i] = uf_hashes_vectors[i].size();
        if (hash_vector_size[i] > 0)
            pq.emplace(make_tuple(uf_hashes_vectors[i][hash_vector_idx[i]++], i));
    }

    uint64_t prev = 0;
    while (pq.size() > 0)
    {
        std::tuple<uint64_t, int> elt = pq.top(); pq.pop();
        uint64_t cur = get<0>(elt);
        //std::cout << "got " << cur << " queue " << get<1>(elt) << std::endl;
        if (cur != prev)
        {
            currentbag->insert(cur);
            nb_elts_pass ++;
        }
        prev = cur;
        
        int i = get<1>(elt);
        if (hash_vector_idx[i] < hash_vector_size[i])
            pq.emplace(make_tuple(uf_hashes_vectors[i][hash_vector_idx[i]++], i));
    }

    for (int i = 0; i < nb_threads; i++)
    free_memory_vector(uf_hashes_vectors[i]);

    
    currentbag->flush();
    
    free_memory_vector(uf_hashes_vectors);

    logging("pass " + to_string(pass+1) + "/" + to_string(nb_passes) + ", " + std::to_string(nb_elts_pass) + " unique hashes written to disk, size " + to_string(nb_elts_pass* sizeof(partition_t) / 1024/1024) + " MB");

    nb_elts += nb_elts_pass;
}


/* main */
template<size_t SPAN>
void bglue(Storage *storage, 
        std::string prefix,
        int kmerSize, 
        int nb_glue_partitions, 
        int nb_threads, 
        bool verbose
        )
{
    auto start_t=chrono::system_clock::now();
    double unit = 1000000000;
    cout.setf(ios_base::fixed);
    cout.precision(1);

    std::cout << "bglue_algo params, prefix:" << prefix << " k:" << kmerSize << " threads:" << nb_threads << std::endl;
    bcalm_logging = verbose;
    size_t k = kmerSize;
    bool debug_uf_stats = false; // formerly cmdline parameter
    bool only_uf = false; // idem

    logging("Starting bglue with " + std::to_string( nb_threads) + " threads");

    //int nbGluePartitions=200; // no longer fixed 
    // autodetecting number of partitions
    int max_open_files = System::file().getMaxFilesNumber() / 2;
    int nbGluePartitions = std::min(2000, max_open_files); // ceil it at 2000 anyhow

    if (nb_glue_partitions > 0)
    {
        nbGluePartitions = nb_glue_partitions;
        logging("Using user-defined number of glue partitions: " + std::to_string( nb_glue_partitions));
    }


    // create a hasher for UF
    typedef typename Kmer<SPAN>::ModelCanonical ModelCanon;
    ModelCanon modelCanon(kmerSize); // i'm a bit lost with those models.. I think GATB could be made more simple here.
    Hasher_T<ModelCanon> hasher(modelCanon);

    ifstream f((prefix + ".glue").c_str());
    if (f.peek() == std::ifstream::traits_type::eof())
    {
        std::cout << "Empty glue file (no sequences)." << std::endl;
        return;
    }

    IBank *in = Bank::open (prefix + ".glue");
    LOCAL(in);
    
    uint64_t nb_glue_sequences = 0;
    
    if (storage != nullptr)
    {
        Group& bcalmGroup = storage->getGroup("bcalm"); 
        nb_glue_sequences = atol(bcalmGroup.getProperty ("nb_sequences_in_glue").c_str());
    }

    if (nb_glue_sequences == 0)
    {
        uint64_t estimated_nb_glue_sequences = in->estimateNbItems();
        logging("estimating number of sequences to be glued (couldn't find true number)");
        nb_glue_sequences = estimated_nb_glue_sequences;
    }
    logging("number of sequences to be glued: "  + to_string(nb_glue_sequences) );

    /*
     * puts all the uf hashes in disk.
     */
    int nb_prepare_passes = 3;
    uint64_t nb_elts = 0;
    for (int pass = 0; pass < nb_prepare_passes; pass++)
        prepare_uf<SPAN>(prefix, in, nb_threads, kmerSize, pass, nb_prepare_passes, nb_elts, nb_glue_sequences);

    // load uf hashes from disk
    std::vector<partition_t> uf_hashes;
    uf_hashes.reserve(nb_elts);
    for (int pass = 0; pass < nb_prepare_passes; pass++)
    {
        IteratorFile<uint64_t> file(prefix+".glue.hashes." + to_string(pass));
        for (file.first(); !file.isDone(); file.next())
            uf_hashes.push_back(file.item());
    }
    if (uf_hashes.size() == 0) // prevent an edge case when there's nothing to glue, boophf doesn't like it
        uf_hashes.push_back(0);

    unsigned long nb_uf_keys = uf_hashes.size();
    logging("loaded all unique UF elements (" + std::to_string(nb_uf_keys) + ") into a single vector of size " + to_string(nb_uf_keys* sizeof(partition_t) / 1024/1024) + " MB");

    int gamma = 3; // make it even faster.
 
    // TODO: why not iterate the UF hashes from disk instead of loading them in memory?
    boomphf::mphf<partition_t , /*TODO we don't need hasher_t here now that we're not hashing kmers, but I forgot to change*/ hasher_t< partition_t> > uf_mphf(nb_uf_keys, uf_hashes, nb_threads, gamma, verbose);

    free_memory_vector(uf_hashes);

    if (verbose)
    {
        unsigned long uf_mphf_memory = uf_mphf.totalBitSize();
        logging("UF MPHF constructed (" + std::to_string(uf_mphf_memory/8/1024/1024) + " MB)" );
    }

    // create a UF data structure
    // this one stores nb_uf_keys * uint64_t (actually, atomic's). so it's bigger than uf_hashes
    unionFind ufkmers(nb_uf_keys);

#if 0
    unionFind<unsigned int> ufmin;
    unionFind<std::string> ufprefixes;
    unsigned int prefix_length = 10;
    unionFind<std::string> ufkmerstr;
#endif
// those were toy one, here is the real one:
    
    // instead of UF of kmers, we do a union find of hashes of kmers. less memory. will have collisions, but that's okay i think. let's see.
    // actually, in the current implementation, partition_t is not used, but values are indeed hardcoded in 32 bits (the UF implementation uses a 64 bits hash table for internal stuff)

    // We loop over sequences.

    /* // uncomment for non-dispatcher version
    auto it = in->iterator();
    for (it->first(); !it->isDone(); it->next())
    {
        const string seq = (*it)->toString();
        const string comment = (*it)->getComment();
    */
    
    auto createUF = [k, &modelCanon, \
        &uf_mphf, &ufkmers, &hasher](const Sequence& sequence)
    {
        const string seq = sequence.toString();
        const string comment = sequence.getComment();

        if (seq.size() < k)
        {
            std::cout << "unexpectedly small sequence found ("<<seq.size()<<"). did you set k correctly?" <<std::endl; exit(1);
        }

        bool lmark = comment[0] == '1';
        bool rmark = comment[1] == '1';

        if ((!lmark) || (!rmark)) // if either mark is 0, no need to associate kmers in UF
            return;

        const string kmerBegin = seq.substr(0, k );
        const string kmerEnd = seq.substr(seq.size() - k , k );

        // UF of canonical kmers in ModelCanon form, then hashed
        const typename ModelCanon::Kmer kmmerBegin = modelCanon.codeSeed(kmerBegin.c_str(), Data::ASCII);
        const typename ModelCanon::Kmer kmmerEnd = modelCanon.codeSeed(kmerEnd.c_str(), Data::ASCII);

        uint32_t v1 = uf_mphf.lookup(hasher(kmmerBegin));
        uint32_t v2 = uf_mphf.lookup(hasher(kmmerEnd));

        ufkmers.union_(v1,v2);
        //ufkmers.union_((hasher(kmmerBegin)), (hasher(kmmerEnd)));

#if 0

        Model::Kmer kmmerBegin = model.codeSeed(kmerBegin.c_str(), Data::ASCII);
        Model::Kmer kmmerEnd = model.codeSeed(kmerEnd.c_str(), Data::ASCII);

        // UF of canonical kmers in string form, not hashed
        string canonicalKmerBegin = modelK1.toString(kmmerBegin.value());
        string canonicalKmerEnd = modelK1.toString(kmmerEnd.value());
        ufkmerstr.union_(canonicalKmerBegin, canonicalKmerEnd);

        // UF of minimizers of kmers
        size_t leftMin(modelK1.getMinimizerValue(kmmerBegin.value()));
        size_t rightMin(modelK1.getMinimizerValue(kmmerEnd.value()));
        ufmin.union_(leftMin, rightMin);

        // UF of prefix of kmers in string form
        string prefixCanonicalKmerBegin = canonicalKmerBegin.substr(0, prefix_length);
        string prefixCanonicalKmerEnd = canonicalKmerEnd.substr(0, prefix_length);
        ufprefixes.union_(prefixCanonicalKmerBegin, prefixCanonicalKmerEnd);
#endif


    };

    Dispatcher dispatcher (nb_threads);
    dispatcher.iterate (in->iterator(), createUF);
    
#if 0
    ufmin.printStats("uf minimizers");

    ufprefixes.printStats("uf " + to_string(prefix_length) + "-prefixes of kmers");

    ufkmerstr.printStats("uf kmers, std::string");
#endif


    logging("UF constructed");

    // remove the UF hashes, could have done that earlier (right before bbhash), but i'd like to keep them in case of segfault
    for (int pass = 0; pass < nb_prepare_passes; pass++)
        System::file().remove (prefix+".glue.hashes." + to_string(pass));


    if (debug_uf_stats) // for debugging
    {
        ufkmers.printStats("uf kmers");
        //ufkmers.dump("uf.dump");
        logging("after computing UF stats");
    }

    if (only_uf) // for debugging
        return;

    /* now we're mirroring the UF to a vector of uint32_t's, it will take less space, and strictly same information
     * this is to get rid of the rank (one uint32) per element in the current UF implementation 
     * we're using the disk to save space of populating one vector from the other in memory. */
    
    BagFile<uint64_t> *ufkmers_bagf = new BagFile<uint64_t>(prefix+".glue.uf");  LOCAL(ufkmers_bagf);
	BagCache<uint64_t> *ufkmers_bag = new BagCache<uint64_t>(  ufkmers_bagf, 10000 );   LOCAL(ufkmers_bag);

    for (unsigned long i = 0; i < nb_uf_keys; i++)
        //ufkmers_vector[i] = ufkmers.find(i); // just in-memory without the disk
        ufkmers_bag->insert(ufkmers.find(i));

    uint64_t size_mdata = sizeof(std::atomic<uint64_t>) * ufkmers.mData.size();
    free_memory_vector(ufkmers.mData);

    logging("freed original UF (" + to_string(size_mdata/1024/1024) + " MB)");

    ufkmers_bag->flush();

    std::vector<uint32_t > ufkmers_vector(nb_uf_keys);
    IteratorFile<uint64_t> ufkmers_file(prefix+".glue.uf");
    unsigned long i = 0;
    for (ufkmers_file.first(); !ufkmers_file.isDone(); ufkmers_file.next())
            ufkmers_vector[i++] = ufkmers_file.item();

    System::file().remove (prefix+".glue.uf");
    
    logging("loaded 32-bit UF (" + to_string(nb_uf_keys*sizeof(uint32_t)/1024/1024) + " MB)");
  
    // setup output file
    string output_prefix = prefix;
    std::atomic<unsigned long> out_id; // identified for output sequences
    out_id = 0;
    BufferedFasta out (output_prefix, 100000);

    auto get_UFclass = [&modelCanon, &ufkmers_vector, &hasher, &uf_mphf]
        (const string &kmerBegin, const string &kmerEnd,
         bool lmark, bool rmark,
         typename ModelCanon::Kmer &kmmerBegin, typename ModelCanon::Kmer &kmmerEnd,  // those will be populated based on lmark and rmark
         bool &found_class)
        {
            found_class = false;
            uint32_t ufclass = 0;

            if (lmark)
            {
                kmmerBegin = modelCanon.codeSeed(kmerBegin.c_str(), Data::ASCII);
                found_class = true;
                ufclass = ufkmers_vector[uf_mphf.lookup(hasher(kmmerBegin))];
            }

            if (rmark)
            {
                kmmerEnd = modelCanon.codeSeed(kmerEnd.c_str(), Data::ASCII);

                if (found_class) // just do a small check
                {
                    if (ufkmers_vector[uf_mphf.lookup(hasher(kmmerEnd))] != ufclass)
                    { std::cout << "bad UF! left kmer has partition " << ufclass << " but right kmer has partition " << ufkmers_vector[uf_mphf.lookup(hasher(kmmerEnd))] << std::endl; exit(1); }
                }
                else
                {
                    ufclass = ufkmers_vector[uf_mphf.lookup(hasher(kmmerEnd))];
                    found_class = true;
                }
            }

            return ufclass;
        };

    std::mutex outLock; // for the main output file
    std::vector<BufferedFasta*> gluePartitions(nbGluePartitions);
    std::string gluePartition_prefix = output_prefix + ".gluePartition.";
    unsigned int max_buffer = 50000;
    std::vector<std::atomic<unsigned long>> nb_seqs_in_partition(nbGluePartitions);


    for (int i = 0; i < nbGluePartitions; i++)
    {
        string filename = gluePartition_prefix + std::to_string(i);
        if (System::file().doesExist(filename))
           System::file().remove (filename);
        gluePartitions[i] = new BufferedFasta(filename, max_buffer);
        nb_seqs_in_partition[i] = 0;
    }

    logging( "Allowed " + to_string((max_buffer * nbGluePartitions) /1024 /1024) + " MB memory for buffers");


    // partition the glue into many files, à la dsk
    std::mutex mtx; // lock to avoid a nasty bug when calling output()
    auto partitionGlue = [k, &modelCanon /* crashes if copied!*/, \
        &get_UFclass, &gluePartitions,
        &out, &outLock, &nb_seqs_in_partition, &out_id, nbGluePartitions, &mtx]
            (const Sequence& sequence)
    {
        const string &seq = sequence.toString();
        const string &comment = sequence.getComment();

        bool lmark = comment[0] == '1';
        bool rmark = comment[1] == '1';

        const string kmerBegin = seq.substr(0, k );
        const string kmerEnd = seq.substr(seq.size() - k , k );

        // make canonical kmer
        typename ModelCanon::Kmer kmmerBegin;
        typename ModelCanon::Kmer kmmerEnd;

        bool found_class = false;

        uint32_t ufclass = get_UFclass(kmerBegin, kmerEnd, lmark, rmark, kmmerBegin, kmmerEnd, found_class);

        if (!found_class) // this one doesn't need to be glued
        {
            const string abundances = comment.substr(3);
            float mean_abundance = get_mean_abundance(abundances);
            uint32_t sum_abundances = get_sum_abundance(abundances);
            
            // for some reason i do need that lock_guard here.. even though output is itself lock guarded. maybe some lazyness in the evauation of the to_string(out_id++)? who kon
            // anyway this fixes the problem, i'll understand it some other time.
            std::lock_guard<std::mutex> lock(mtx);
            output(seq, out, std::to_string(out_id++) + " LN:i:" + to_string(seq.size()) + " KC:i:" + to_string(sum_abundances) + " km:f:" + to_string_with_precision(mean_abundance)); 
            // km is not a standard GFA field so i'm putting it in lower case as per the spec
            // maybe could optimize by writing to disk using queues, if that's ever a bottleneck
            return;
        }

        int index = ufclass % nbGluePartitions;
        //stringstream ss1; // to save partition later in the comment. [why? probably to avoid recomputing it]
        //ss1 << blabla;

        output(seq, *gluePartitions[index], comment);
        nb_seqs_in_partition[index]++;
    };

    logging("Disk partitioning of glue");
    dispatcher.iterate (in->iterator(), partitionGlue); // multi-threaded
    /*// single-threaded version
     auto it = in->iterator();    
     for (it->first (); !it->isDone(); it->next())
        partitionGlue(it->item());
    */

    for (int i = 0; i < nbGluePartitions; i++)
        delete gluePartitions[i]; // takes care of the final flush (this doesn't delete the file, just closes it)
    free_memory_vector(gluePartitions);
    out.flush();
 

    logging("Done disk partitioning of glue");

    // get top10 largest glue partitions
    int top_n_glue_partition = std::min(10,nbGluePartitions);
    vector<unsigned long> vx, copy_nb_seqs_in_partition;
    vx.resize(nb_seqs_in_partition.size());
    copy_nb_seqs_in_partition.resize(nb_seqs_in_partition.size());
    for(unsigned int i= 0; i<nb_seqs_in_partition.size(); ++i ) 
    {
        vx[i]= i;
        copy_nb_seqs_in_partition[i] = nb_seqs_in_partition[i]; // to get rid of atomic type
    }
    partial_sort( vx.begin(), vx.begin()+top_n_glue_partition, vx.end(), Comp<unsigned long>(copy_nb_seqs_in_partition) );

    if (verbose)
    {
        std::cout << "Top 10 glue partitions by size:" << std::endl;
        for (int i = 0; i < top_n_glue_partition; i++)
            std::cout << "Glue partition " << vx[i] << " has " << copy_nb_seqs_in_partition[vx[i]] << " sequences " << endl;
    }

    logging("Glueing partitions");

    // glue all partitions using a thread pool
    ThreadPool pool(nb_threads);
    for (int partition = 0; partition < nbGluePartitions; partition++)
    {
        auto glue_partition = [&modelCanon, &ufkmers, partition, &gluePartition_prefix, nbGluePartitions, &copy_nb_seqs_in_partition,
        &get_UFclass, &out, &outLock, &out_id, kmerSize, &mtx]( int thread_id)
        {
            int k = kmerSize;

            string partitionFile = gluePartition_prefix + std::to_string(partition);
            BankFasta partitionBank (partitionFile); // BankFasta
            BankFasta::Iterator it (partitionBank); // BankFasta

            outLock.lock(); // should use a printlock..
            if (partition % 20 == 0) // sparse printing
            {
                string message = "Gluing partition " +to_string(partition) + " (size: " +to_string(System::file().getSize(partitionFile)/1024/1024) + " MB)";
                logging(message);
            }
            outLock.unlock();

            unordered_map<int, vector< markedSeq<SPAN> >> msInPart;
            uint64_t seq_index = 0;

            for (it.first(); !it.isDone(); it.next()) // BankFasta
            {
                const string seq = it->toString();
                const string comment = it->getComment();

                const string kmerBegin = seq.substr(0, k );
                const string kmerEnd = seq.substr(seq.size() - k , k );

                uint32_t ufclass = 0;
                bool found_class = false;

                bool lmark = comment[0] == '1';
                bool rmark = comment[1] == '1';

                // todo speed improvement: get partition id from sequence header (so, save it previously)

                // make canonical kmer
                typename ModelCanon::Kmer kmmerBegin, kmmerEnd;

                ufclass = get_UFclass(kmerBegin, kmerEnd, lmark, rmark, kmmerBegin, kmmerEnd, found_class);

                // compute kmer extremities if we have not already
                if (!lmark)
                    kmmerBegin = modelCanon.codeSeed(kmerBegin.c_str(), Data::ASCII);
                if (!rmark)
                    kmmerEnd = modelCanon.codeSeed(kmerEnd.c_str(), Data::ASCII);

                markedSeq<SPAN> ms(seq_index, lmark, rmark, kmmerBegin.value(), kmmerEnd.value());

                //if (ufclass == 38145) std::cout << " ufclass " << ufclass << " seq " << seq << " seq index " << seq_index << " " << lmark << rmark << " ks " << kmmerBegin.value() << " ke " << kmmerEnd.value() << std::endl; // debug specific partition
                msInPart[ufclass].push_back(ms);
                seq_index++;
            }

            // now iterates all sequences in a partition to determine the order in which they're going to be glues (avoid intermediate gluing)
            vector<vector<uint32_t>> ordered_sequences_idxs ;
            for (auto it = msInPart.begin(); it != msInPart.end(); it++)
            {
                bool debug = false; //debug = it->first == 38145; // debug specific partition
                //std::cout << "1.processing partition " << it->first << std::endl;
                determine_order_sequences<SPAN>(ordered_sequences_idxs, it->second, kmerSize, debug); // return indices of markedSeq's inside it->second
                //std::cout << "2.processing partition " << it->first << " nb ordered sequences: " << ordered_sequences_idxs.size() << std::endl;
                free_memory_vector(it->second);
            }

            msInPart.clear();
            unordered_map<int,vector<markedSeq<SPAN>>>().swap(msInPart); // free msInPart
            
            vector<string> sequences;
            vector<string> abundances;
            sequences.reserve(copy_nb_seqs_in_partition[partition]);
            abundances.reserve(copy_nb_seqs_in_partition[partition]);
            
            for (it.first(); !it.isDone(); it.next()) // BankFasta
            {
                const string seq = it->toString();
                const string comment = it->getComment();
                const string abundance_str = comment.substr(3);
                sequences.push_back(seq);
                abundances.push_back(abundance_str);
            }

            for (auto itO = ordered_sequences_idxs.begin(); itO != ordered_sequences_idxs.end(); itO++)
            {
                string seq, abs;
                glue_sequences(*itO, sequences, abundances, kmerSize, seq, abs); // takes as input the indices of ordered sequences, and the markedSeq's themselves

                float mean_abundance = get_mean_abundance(abs);
                uint32_t sum_abundances = get_sum_abundance(abs);
                {
                    // for some reason i do need that lock_guard here.. even though output is itself lock guarded. maybe some lazyness in the evauation of the to_string(out_id++)? who kon
                    // anyway this fixes the problem, i'll understand it some other time.
                    std::lock_guard<std::mutex> lock(mtx);
                    output(seq, out, std::to_string(out_id++) + " LN:i:" + to_string(seq.size()) + " KC:i:" + to_string(sum_abundances) + " km:f:" + to_string_with_precision(mean_abundance)); 
                }
            }
                
            free_memory_vector(ordered_sequences_idxs);

            partitionBank.finalize(); // BankFasta

            System::file().remove (partitionFile);

        };

        pool.enqueue(glue_partition);
        //glue_partition(0); // single threaded
    }

    pool.join();
   
   out.flush(); // not sure if necessary

    logging("end");

    bool debug_keep_glue_files = true; // for debugging // TODO enable it if -redo-bglue param was provided (need some info from UnitigsConstructionAlgorithm). 
    if (debug_keep_glue_files)
    {
        std::cout << "debug: not deleting glue files" << std::endl;
    }
    else
    {
        // cleanup glue files
        std::string line;
        std::ifstream infile(prefix + ".glue");
        while (std::getline(infile, line))
        {
            System::file().remove (line); 
        }
        infile.close();
        System::file().remove (prefix + ".glue");
    }
    auto end_t=chrono::system_clock::now();
    float wtime = chrono::duration_cast<chrono::nanoseconds>(end_t - start_t).count() / unit;

    if (storage != nullptr)
    {
        Group& bcalmGroup = storage->getGroup("bcalm"); 
        bcalmGroup.setProperty ("wtime_glue",     Stringify::format("%f", wtime));
    }
}

}}}}
