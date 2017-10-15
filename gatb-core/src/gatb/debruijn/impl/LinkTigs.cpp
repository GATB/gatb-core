#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bcalm2/logging.hpp>
#include <gatb/debruijn/impl/ExtremityInfo.hpp>
#include <gatb/debruijn/impl/LinkTigs.hpp>
#include <gatb/kmer/impl/Model.hpp> // for revcomp_4NT

#include <queue>
#include <string>
#include <unordered_map>


using namespace std;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::debruijn;
using namespace gatb::core::debruijn::impl;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

namespace gatb { namespace core { namespace debruijn { namespace impl  {

    static constexpr int nb_passes = 8;
    static void write_final_output(const string& unitigs_filename, bool verbose, BankFasta* out, uint64_t &nb_unitigs);
    static bool get_link_from_file(std::ifstream& input, std::string &link, uint64_t &unitig_id);

/* this procedure finds the overlaps between unitigs, using a hash table of all extremity (k-1)-mers
 * I guess it's like AdjList in ABySS. It's also like contigs_to_fastg in MEGAHIT.
 * 
 * could be optimized by keeping edges during the BCALM step and tracking kmers in unitigs, but it's not the case for now, because would need to modify ograph 
 *
 * it uses the disk to store the links for extremities until they're merged into the final unitigs file. 
 *
 * so the memory usage is just that of the hash tables that record kmers, not of the links
 */
template<size_t span>
void link_tigs(string unitigs_filename, int kmerSize, int nb_threads, uint64_t &nb_unitigs, bool verbose)
{
    bcalm_logging = verbose;
    BankFasta* out = new BankFasta(unitigs_filename+".indexed");
    if (kmerSize < 4) { std::cout << "error, recent optimizations (specifically link_unitigs) don't support k<5 for now" << std::endl; exit(1); }
    logging("Finding links between unitigs");

    for (int pass = 0; pass < nb_passes; pass++)
        link_unitigs_pass<span>(unitigs_filename, verbose, pass, kmerSize);

    write_final_output(unitigs_filename, verbose, out, nb_unitigs);
   
    delete out;
    system::impl::System::file().remove (unitigs_filename);
    system::impl::System::file().rename (unitigs_filename+".indexed", unitigs_filename);

    logging("Done finding links between unitigs");
}


/*
 * takes all the prefix.links.* files, sorted by unitigs.
 * do a n-way merge to gather the links for each unitig in unitig order
 * (single-threaded)
 */


static void write_final_output(const string& unitigs_filename, bool verbose, BankFasta* out, uint64_t &nb_unitigs)
{
    logging("gathering links from disk");
    std::ifstream* inputLinks[nb_passes];

    bitset<nb_passes> finished;
    typedef std::tuple<uint64_t /*unitig id*/, int /*pass*/, std::string /*links */> pq_elt_t;
    priority_queue<pq_elt_t, vector<pq_elt_t>, std::greater<pq_elt_t> > pq;
    
    BankFasta inputBank (unitigs_filename);
    BankFasta::Iterator itSeq (inputBank);
    itSeq.first(); 
    string cur_links, seq, comment;
    seq = itSeq->toString();
    comment = itSeq->getComment();
 
    for (int pass = 0; pass < nb_passes; pass++)
    {
        string link; uint64_t unitig;
        inputLinks[pass] = new std::ifstream(unitigs_filename+ ".links." + to_string(pass));
        // prime the pq with the first element in the file
        if (get_link_from_file(*inputLinks[pass], link, unitig))
            pq.emplace(make_tuple(unitig, pass, link));
        else
            finished[pass] = true;
    }

    uint64_t last_unitig = 0;
    nb_unitigs = 0; // passed variable

    // nb_passes-way merge sort
    while ((!finished.all()) || pq.size() > 0)
    {
        pq_elt_t cur = pq.top(); pq.pop();
        int pass = get<1>(cur);
        uint64_t unitig = get<0>(cur);

        if (unitig != last_unitig)
        {
            Sequence s (Data::ASCII);
            s.getData().setRef ((char*)seq.c_str(), seq.size());
            s._comment = comment + " " + cur_links;
            out->insert(s);
            
            cur_links = "";
            nb_unitigs++;
            last_unitig = unitig;
            itSeq.next();
            seq = itSeq->toString();
            comment = itSeq->getComment();
        }
            
        cur_links += get<2>(cur);
        //if (unitig < 10)  std::cout << " popped " << pass << " " << unitig << " " << cur_links << std::endl; // debug

        // read next entry in the inputLinks[pass] file that we just popped
        if (finished[pass])
            continue;
        string link;
        if (get_link_from_file(*inputLinks[pass], link, unitig))
            pq.emplace(make_tuple(unitig, pass, link));
        else
            finished[pass] = true;

    }
    // write the last element
    Sequence s (Data::ASCII);
    s.getData().setRef ((char*)seq.c_str(), seq.size());
    s._comment = comment + " " + cur_links;
    out->insert(s);
    nb_unitigs++;

    for (int pass = 0; pass < nb_passes; pass++)
    {
        system::impl::System::file().remove (unitigs_filename + ".links." + to_string(pass));
        delete inputLinks[pass];
    }
}


// well well, some potential code duplication with Model.hpp in here (or rather, specialization), but sshh
static inline int nt2int(char nt)
{
    if (nt=='A') return 0;
    if (nt=='C') return 1;
    if (nt=='T') return 2;
    if (nt=='G') return 3;
    return 0;
}

/* that code doesn't support more than 8 passes*/
static int normalized_smallmer(const unsigned char c1, const unsigned char c2, const unsigned char c3, const unsigned char c4)
{
    unsigned char smallmer = (nt2int(c1)<<6) + (nt2int(c2)<<4) + (nt2int(c3)<<2) + nt2int(c4);
    const unsigned char rev = revcomp_4NT[smallmer];
    if (rev < smallmer)
        smallmer = rev;
    return smallmer;
}

static bool
is_in_pass (const std::string &seq, int pass, Unitig_pos p, int kmerSize) // TODO this is so un-even. should do more proper hashing..
{
    int e = 0;
    if (p == UNITIG_END)
        e = seq.size()-(kmerSize-1);
    // x = 0123456789
    // k = 5, k-1=4
    // seq.size()-1-(k-1) = 10-4 = 6
    return (normalized_smallmer(seq[e],seq[e+1],seq[e+kmerSize-1-1-1],seq[e+kmerSize-1-1]) % nb_passes) == pass;
}

/* returns true if it has read an element, false otherwise */
static bool get_link_from_file(std::ifstream& input, std::string &link, uint64_t &unitig_id)
{
    string line;
    if (std::getline(input, line))
    {
        unitig_id = stoull(line);
    }
    else
        return false;
    if (std::getline(input, link))
    {
    }
    else
        return false;
    return true;
}


static void record_links(uint64_t utig_id, int pass, const string &link, std::ofstream &links_file)
{
    // maybe do a buffered thing here but it's not clear if it is bottleneck. in CAMI-medium it took 6 mins without if i don't write the links_file and 8 mins if i do..
    links_file << to_string(utig_id) << "\n";
    links_file << link << "\n";
}


template<size_t span>
void link_unitigs_pass(const string unitigs_filename, bool verbose, const int pass, const int kmerSize)
{
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::Type           Type;

    bool debug = false;

    BankFasta inputBank (unitigs_filename);
    BankFasta::Iterator itSeq (inputBank);
    uint64_t utig_counter = 0;
    
    Model modelKminusOne(kmerSize - 1); // it's canonical (defined in the .hpp file)
    
    typedef typename std::unordered_map<Type, std::vector<uint64_t>> NodeLinksMap;

    NodeLinksMap utigs_links_map;

    logging("step 1 pass " + to_string(pass));

    // this is the memory-limiting step, but can be lowered with larger nb_pass
    for (itSeq.first(); !itSeq.isDone(); itSeq.next()) 
    {
        const string& seq = itSeq->toString();
        if (debug) std::cout << "unitig: " << seq << std::endl;

        if (is_in_pass(seq, pass, UNITIG_BEGIN, kmerSize))
        { 
            if (debug) std::cout << "pass " << pass << " examining beginning" << std::endl;
            typename Model::Kmer kmerBegin = modelKminusOne.codeSeed(seq.substr(0, kmerSize-1).c_str(), Data::ASCII);
            bool beginInSameOrientation = modelKminusOne.toString(kmerBegin.value()) == seq.substr(0,kmerSize-1);
            ExtremityInfo eBegin(utig_counter, !beginInSameOrientation /* because we record rc*/, UNITIG_BEGIN);
            utigs_links_map[kmerBegin.value()].push_back(eBegin.pack());
        }
        if (is_in_pass(seq, pass, UNITIG_END, kmerSize))
        {
            if (debug) std::cout << "pass " << pass << " examining end" << std::endl;
            typename Model::Kmer kmerEnd = modelKminusOne.codeSeed(seq.substr(seq.size() - kmerSize+1).c_str(), Data::ASCII);
            bool endInSameOrientation = modelKminusOne.toString(kmerEnd.value()) == seq.substr(seq.size() - kmerSize+1);
            ExtremityInfo eEnd(  utig_counter, !endInSameOrientation,                             UNITIG_END);
            utigs_links_map[kmerEnd.value()].push_back(eEnd.pack());
            // there is no UNITIG_BOTH here because we're taking (k-1)-mers.
        }
        utig_counter++;
    }

    std::ofstream links_file(unitigs_filename+".links." +to_string(pass));

    uint64_t nb_hashed_entries = 0;
    for (auto v : utigs_links_map)
        nb_hashed_entries += v.second.size(); 
    logging("step 2 (" + to_string(utigs_links_map.size()) + "kmers/" + to_string(nb_hashed_entries) + "extremities)");

    utig_counter = 0;
    for (itSeq.first(); !itSeq.isDone(); itSeq.next()) 
    {
        const string& seq = itSeq->toString();
        
        if (debug) std::cout << "unitig: " << seq << std::endl;
 
        if (is_in_pass(seq, pass, UNITIG_BEGIN, kmerSize))
        {
            if (debug) std::cout << "pass " << pass << " examining beginning" << std::endl;
            typename Model::Kmer kmerBegin = modelKminusOne.codeSeed(seq.substr(0, kmerSize-1).c_str(), Data::ASCII);
            bool beginInSameOrientation =  modelKminusOne.toString(kmerBegin.value()) == seq.substr(0,kmerSize-1); // that could be optimized, revcomp was already computed during codeSeed
            // treat special palindromic kmer cases
            bool nevermindInOrientation = false;
            if (((kmerSize - 1) % 2) == 0) 
                if (kmerBegin.isPalindrome()) nevermindInOrientation = true;

            string in_links = " "; // necessary placeholder to indicate we have links for that unitig

            // in-neighbors
            // so, the current kmer is kmerBegin, and we're iterating all potential in-neighbors as e_in objects
            for (auto in_packed : utigs_links_map[kmerBegin.value()])
            {
                ExtremityInfo e_in(in_packed);

                if (debug) std::cout << "extremity " << modelKminusOne.toString(kmerBegin.value()) << " ";
                if (debug) std::cout << "potential in-neighbor: " << e_in.toString() << " beginSameOrientation " << beginInSameOrientation;

                // what we want are these four cases:
                //  ------[end same orientation] -> [begin same orientation]----
                //  [begin diff orientation]---- -> [begin same orientation]----
                //  ------[end diff orientation] -> [begin diff orientation]----
                //  [begin same orientation]---- -> [begin diff orientation]----
                if ((       ((beginInSameOrientation) &&    (e_in.pos == UNITIG_END  ) && (e_in.rc == false)) ||
                            ((beginInSameOrientation) &&    (e_in.pos == UNITIG_BEGIN) && (e_in.rc == true)) ||
                            (((!beginInSameOrientation)) && (e_in.pos == UNITIG_END  ) && (e_in.rc == true)) ||
                            (((!beginInSameOrientation)) && (e_in.pos == UNITIG_BEGIN) && (e_in.rc == false)))
                        || nevermindInOrientation)
                {
                    
                    // this was for when i was wanting to save space while storing links in memory. now storing on disk
                    //LinkInfo li(e_in.unitig, e_in.rc ^ beginInSameOrientation);
                    //incoming[utig_number].push_back(li.pack());

                    //bool rc = e_in.rc ^ (!beginInSameOrientation); // "rc" sets the destination strand // i don't think it's the right formula because of k-1-mers that are their self revcomp. see the mikko bug in the test folder, that provides a nice illustration of that 
                    bool rc = e_in.pos == UNITIG_END; // a better way to determine the rc flag is just looking at position of e_in k-1-mer 

                    in_links += "L:-:" + to_string(e_in.unitig) + ":" + (rc?"-":"+") + " "; 

                    /* what to do when kmerBegin is same as forward and reverse?
                     used to have this:
                        if (nevermindInOrientation)
                          in_links += "L:-:" + to_string(e_in.unitig) + ":" + ((!rc)?"+":"-") + " "; // in that case, there is also another link with the reverse direction
                     but it's bogus, because the reverse direction of the other sequence won't have the same start k-1-mer, even if it is just k-long
                     so there isn't anything to do actually. the case is handled by the "|| nevermindInOrientation )" in the if above */
                     

                    if (debug) std::cout << " [valid] ";
                }
                if (debug) std::cout << std::endl;
            }
            
            record_links(utig_counter, pass, in_links, links_file);
        }


        if (is_in_pass(seq, pass, UNITIG_END, kmerSize))
        {
            if (debug) std::cout << "pass " << pass << " examining end" << std::endl;
            typename Model::Kmer kmerEnd = modelKminusOne.codeSeed(seq.substr(seq.size() - kmerSize+1).c_str(), Data::ASCII);
            bool endInSameOrientation =  modelKminusOne.toString(kmerEnd.value()) == seq.substr(seq.size() - kmerSize+1);

            bool nevermindOutOrientation = false;
            if (((kmerSize - 1) % 2) == 0) 
                if (kmerEnd.isPalindrome())   nevermindOutOrientation = true;

            string out_links = " "; // necessary placeholder to indicate we have links for that unitig

            // out-neighbors
            for (auto out_packed : utigs_links_map[kmerEnd.value()])
            {
                ExtremityInfo e_out(out_packed);

                if (debug) std::cout << "extremity " << modelKminusOne.toString(kmerEnd.value()) << " ";
                if (debug) std::cout << "potential out-neighbor: " << e_out.toString();

                // what we want are these four cases:
                //  ------[end same orientation] -> [begin same orientation]----
                //  ------[end same orientation] -> ------[end diff orientation]
                //  ------[end diff orientation] -> [begin diff orientation]----
                //  ------[end diff orientation] -> ------[end same orientation]
                if ((((endInSameOrientation) && (e_out.pos == UNITIG_BEGIN) && (e_out.rc == false)) ||
                            ((endInSameOrientation) && (e_out.pos == UNITIG_END  ) && (e_out.rc == true)) ||
                            (((!endInSameOrientation)) && (e_out.pos == UNITIG_BEGIN) && (e_out.rc == true)) ||
                            (((!endInSameOrientation)) && (e_out.pos == UNITIG_END  ) && (e_out.rc == false)))
                        ||nevermindOutOrientation)
                {

                    //LinkInfo li(e_out.unitig, e_out.rc ^ endInSameOrientation);
                    //outcoming[utig_number].push_back(li.pack());
                    
                    bool rc = e_out.pos == UNITIG_END; // a better way to determine the rc flag is just looking at position of e_in k-1-mer 

                    out_links += "L:+:" + to_string(e_out.unitig) + ":" + (rc?"-":"+") + " "; 

                    if (debug) std::cout << " [valid] ";
                }
                if (debug) std::cout << std::endl;
            }
            record_links(utig_counter, pass, out_links, links_file);
        }

        utig_counter++;
    }
}

}}}}
