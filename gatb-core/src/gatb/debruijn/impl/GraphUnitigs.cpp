/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <gatb/debruijn/impl/GraphUnitigs.hpp>

#include <gatb/bank/api/IBank.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/system/api/IThread.hpp> // for ISynchronizer 

#include <gatb/tools/collections/impl/ContainerSet.hpp>
#include <gatb/tools/collections/impl/IterableHelpers.hpp>

#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/LibraryInfo.hpp>
#include <gatb/tools/misc/impl/HostInfo.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/tools/misc/impl/Tool.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/CountProcessor.hpp>
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>

#include <gatb/debruijn/impl/Simplifications.hpp>

using namespace std;

using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::debruijn::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

#undef NDEBUG
#include <cassert>

#define DEBUG(a)  //a

// some shorthands for unitigs
#define BaseGraph GraphTemplate< NodeFast<span>, EdgeFast<span>, GraphDataVariantFast<span> >

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPHUNITIGSCPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPHUNITIGSCPP_

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

/********************************************************************************
                 #####   ######      #     ######   #     #
                #     #  #     #    # #    #     #  #     #
                #        #     #   #   #   #     #  #     #
                #  ####  ######   #     #  ######   #######
                #     #  #   #    #######  #        #     #
                #     #  #    #   #     #  #        #     #
                 #####   #     #  #     #  #        #     #
********************************************************************************/

// TODO: maybe change a few things in here or else delete and use Graph::getOptionsParser
template<size_t span>
IOptionsParser* GraphUnitigsTemplate<span>::getOptionsParser (bool includeMandatory)
{

    /** We build the root options parser. */
    OptionsParser* parser = new OptionsParser ("graph");

    /** We add children parser to it (kmer count, bloom/debloom, branching). */
    parser->push_back (SortingCountAlgorithm<>::getOptionsParser(includeMandatory));

    /** We create a "general options" parser. */
    IOptionsParser* parserGeneral  = new OptionsParser ("general");
    parserGeneral->push_front (new OptionOneParam (STR_INTEGER_PRECISION, "integers precision (0 for optimized value)", false, "0", false));
    parserGeneral->push_front (new OptionOneParam (STR_VERBOSE,           "verbosity level",      false, "1"  ));
    parserGeneral->push_front (new OptionOneParam (STR_NB_CORES,          "number of cores",      false, "0"  ));
    parserGeneral->push_front (new OptionNoParam  (STR_CONFIG_ONLY,       "dump config only"));

    /** We add it to the root parser. */
    parser->push_back  (parserGeneral);

    return parser;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>  GraphUnitigsTemplate<span>::create (bank::IBank* bank, const char* fmt, ...)
{
    IOptionsParser* parser = getOptionsParser (false);   LOCAL(parser);

    /** We build the command line from the format and the ellipsis. */
    std::string commandLine;
    char* buffer = 0;
    va_list args;
    va_start (args, fmt);
    vasprintf (&buffer, fmt, args);
    va_end (args);
    if (buffer != NULL)  {  commandLine = buffer;  FREE (buffer);  }

    try
    {
        return  GraphUnitigsTemplate(bank, parser->parseString(commandLine));
    }
    catch (OptionFailure& e)
    {
        e.displayErrors (std::cout);

        throw system::Exception ("Graph construction failure because of bad parameters (notify a developer)");
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>  GraphUnitigsTemplate<span>::create (const char* fmt, ...)
{
    IOptionsParser* parser = getOptionsParser (true);   LOCAL (parser);

    /** We build the command line from the format and the ellipsis. */
    std::string commandLine;
    char* buffer = 0;
    va_list args;
    va_start (args, fmt);
    vasprintf (&buffer, fmt, args);
    va_end (args);
    if (buffer != NULL)  {  commandLine = buffer;  FREE (buffer);  }

    try
    {
        return  GraphUnitigsTemplate (parser->parseString(commandLine), true); /* will call the GraphUnitigsTemplate<span>::GraphUnitigsTemplate (tools::misc::IProperties* params, bool load_unitigs_after) constructor */
    }
    catch (OptionFailure& e)
    {
        e.displayErrors (std::cout);
        throw system::Exception ("Graph construction failure because of bad parameters (notify a developer)");
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE : 
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : load a graph from I don't know where. looks like dummy?
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>::GraphUnitigsTemplate (size_t kmerSize)
    : GraphTemplate<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span> >(kmerSize)
{
    // will call Graph's constructor for (kmerSize), no big deal
    //std::cout << "kmersize graphUtemplate constructor" << std::endl;
}

/*********************************************************************
** METHOD  :
** PURPOSE : loads a graph from a h5 file name
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>::GraphUnitigsTemplate (const std::string& uri)
{
    std::cout << "unitigs graph constructor(uri) not supported" << std::endl; exit(1);
}

/*********************************************************************
** METHOD  :
** PURPOSE : creates a graph from a bank
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>::GraphUnitigsTemplate (bank::IBank* bank, tools::misc::IProperties* params)
{
    // quick hack,  not supposed to be used outside o tests
    
    /** We get the kmer size from the user parameters. */
    BaseGraph::_kmerSize = params->getInt (STR_KMER_SIZE);
    modelK = new Model(BaseGraph::_kmerSize);
    modelKdirect= new ModelDirect(BaseGraph::_kmerSize);
    size_t integerPrecision = params->getInt (STR_INTEGER_PRECISION);
    /** We configure the data variant according to the provided kmer size. */
    BaseGraph::setVariant (BaseGraph::_variant, BaseGraph::_kmerSize, integerPrecision);
    string unitigs_filename = "dummy.unitigs.fa"; // because there's already a bank, but we don't know its name maybe? so just to be safe, i'm setting a dummy unitigs file. anyway, this constructor is only called in tests i think, not by minia for sure.

        /*(params->get(STR_URI_OUTPUT) ?
            params->getStr(STR_URI_OUTPUT) :                                                                                                                                                                 System::file().getBaseName (input)
            )+ ".unitigs.fa";   */       

    params->setInt(STR_REPARTITION_TYPE, 1);
    params->setInt(STR_MINIMIZER_TYPE, 1);

    /** We build the graph according to the wanted precision. */
    boost::apply_visitor ( build_visitor_solid<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>>(*this, bank,params),  *(GraphDataVariantFast<span>*)BaseGraph::_variant);

    build_unitigs_postsolid(unitigs_filename, params);

    load_unitigs(unitigs_filename);
}

static /* important that it's static! else TemplateSpecialization8 will complain*/
char revcomp (char s) {
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

static string revcomp (const string &s) {
	string rc;
	for (signed int i = s.length() - 1; i >= 0; i--) {rc += revcomp(s[i]);}
	return rc;
}

template<size_t span>
void GraphUnitigsTemplate<span>::build_unitigs_postsolid(std::string unitigs_filename, tools::misc::IProperties* props)
{
    /** We may have to stop just after configuration. I don't know if that happens in GraphU though. */
    if (props->get(STR_CONFIG_ONLY))  { std::cout << "GraphU Config_only! does that happen?" << std::endl; return; }

    //if (!BaseGraph::checkState(BaseGraph::STATE_SORTING_COUNT_DONE))
    if (!checkState(STATE_SORTING_COUNT_DONE))
    {
        //throw system::Exception ("Graph construction failure during build_visitor_postsolid, the input h5 file needs to contain at least solid kmers.\nIf this is an old .h5 file, created with GATB-Core's Graph instead of GraphUnitigs, please re-create it.");
        // let's try with shared state.
        throw system::Exception ("Graph construction failure during build_visitor_postsolid, the input h5 file needs to contain at least solid kmers.");
    }
    
    bool force_loading_unitigs = false; // debug option

    if (!checkState(STATE_BCALM2_DONE) && (!force_loading_unitigs /* for debug, if unitigs are made but the h5 didn't register it, stupid h5*/))
    {
        int nb_threads =
            props->getInt(STR_NB_CORES);

        size_t  kmerSize = BaseGraph::getKmerSize();
        if (kmerSize != (unsigned int)props->getInt(STR_KMER_SIZE))
            std::cout << "kmer discrepancy: should i take " << kmerSize << " or " << props->getInt(STR_KMER_SIZE) << std::endl;

        UnitigsConstructionAlgorithm<span> unitigs_algo(BaseGraph::getStorage(), unitigs_filename, nb_threads, props);

        BaseGraph::executeAlgorithm(unitigs_algo, &BaseGraph::getStorage(), props, BaseGraph::_info);
    
        nb_unitigs = unitigs_algo.nb_unitigs;
        
        setState(STATE_BCALM2_DONE);
    }
    else
        nb_unitigs = atol (BaseGraph::getGroup().getProperty ("nb_unitigs").c_str());

    /** We save the state at storage root level. */
    BaseGraph::getGroup().setProperty ("state",          Stringify::format("%d", BaseGraph::_state));
    BaseGraph::getGroup().setProperty ("nb_unitigs",     Stringify::format("%d", nb_unitigs));
}

static void
parse_unitig_header(string header, float& mean_abundance, vector<uint64_t>& inc, vector<uint64_t>& outc)
{
    bool debug = false;
    if (debug) std::cout << "parsing unitig links for " << header << std::endl;
    std::stringstream stream(header);
    while(1) {
        string tok;
        stream >> tok;
        if(!stream)
            break;

        if (tok.size() < 3)
            // that's the id, skip it
            continue;

        string field = tok.substr(0,2);
        if (field == "L:")
        {
            bool in = tok.substr(2,1) == "-";
            int pos_rc = tok.find_last_of(':');
            bool rc = tok.substr(pos_rc+1) == "-";
            tok = tok.substr(0,pos_rc); // chop last field
            uint64_t unitig = atoi(tok.substr(tok.find_last_of(':')+1).c_str());
            /* situation is: 
             * L:+:next_unitig:+    unitig[end] -> [begin]next_unitig 
             * L:+:next_unitig:-    unitig[end] -> [begin]next_unitig_rc
             * L:-:next_unitig:+    unitig_rc[end] -> [begin]next_unitig     or alternatively, next_unitig_rc[end] -> [begin]unitig
             * L:-:next_unitig:-    unitig_rc[end] -> [begin]next_unitig_rc                       next_unitig[end] -> [begin]unitig
             * */

            /* setting pos:
             * in case of single-kmer unitig, pos will be wrong (should be UNITIG_BOTH, but i'm not storing this info in just 1 bit). Instead of encoding it here (would add 1 bit), getEdges as well as simplePath_avance will be inferring that it's UNITIG_BOTH in cases where the unitig is just of length k
             * thus, pos is actually also given by the following formula, if you think hard about it and look at the situations above (actually i got super confused and changed this code so many times until all unit tests passed)*/
            Unitig_pos pos = (rc)?UNITIG_END:UNITIG_BEGIN;
            if (in)
                rc = !rc;

            ExtremityInfo li(unitig, rc, pos);
            if (debug)
                std::cout << "inserting "<< (in?"incoming":"outcoming") <<  " extremity " << li.toString() << std::endl;
            if (in)
                inc.push_back(li.pack());
            else
                outc.push_back(li.pack());
        }
        else
        {
            if (field == "KM")
            {
                mean_abundance = atof(tok.substr(tok.find_last_of(':')+1).c_str());
                //std::cout << "unitig " << header << " mean abundance " << mean_abundance << std::endl;
            }
            // we don't care about other fields
        }
    }
}
        

static void
insert_navigational_vector(std::vector<uint64_t> &v, std::vector<uint64_t>& to_insert, uint64_t utig_counter, std::vector<uint64_t> &v_map)
{
    v_map[utig_counter] = v.size();
    v.insert(v.end(), to_insert.begin(), to_insert.end());
}

//http://stackoverflow.com/questions/30540101/iterator-for-a-subset-of-a-vector
template <class Iter>
class range {
    Iter b;
    Iter e;
    public:
    range(Iter b, Iter e) : b(b), e(e) {}
    Iter begin() { return b; }
    Iter end() { return e; }
};

template <class Container>
range<typename Container::const_iterator> 
make_range(Container& c, size_t b, size_t e) {
    return range<typename Container::const_iterator> (c.begin()+b, c.begin()+e);
}

/* returns an iterator of all incoming or outcoming edges from an unitig */
static 
range<std::vector<uint64_t>::const_iterator >
get_from_navigational_vector(const std::vector<uint64_t> &v, uint64_t utig, const std::vector<uint64_t> &v_map) 
{
    if (utig == v_map.size() /*total number of unitigs*/ - 1)
    {
        //std::cout << "get from nav vector " << to_string(utig) << " " << to_string(v_map[utig]) << " " <<  to_string(v.size()) << " last unitig" << std::endl;
        return make_range(v,v_map[utig],v.size());
    }
    else
    {
        //std::cout << "get from nav vector " << to_string(utig) << " " << to_string(v_map[utig]) << " " <<  to_string(v_map[utig+1]) << " (utig " << utig << "/" << v_map.size() << ")" << std::endl;
        return make_range(v,v_map[utig],v_map[utig+1]);
    }
}



template<size_t span>
void GraphUnitigsTemplate<span>::load_unitigs(string unitigs_filename)
{
    BankFasta inputBank (unitigs_filename);
    //bank::IBank* inputBank = Bank::open (unitigs_filename);
    //LOCAL (inputBank);
    //ProgressIterator<bank::Sequence> itSeq (*inputBank, "loading unitigs");
    BankFasta::Iterator itSeq (inputBank);

 
    incoming_map.resize(nb_unitigs);
    outcoming_map.resize(nb_unitigs);

   unsigned int kmerSize = BaseGraph::_kmerSize;

    nb_unitigs_extremities = 0; // will be used by NodeIterator (getNodes)
    uint64_t utig_counter = 0;
    uint64_t nb_utigs_nucl = 0;
    uint64_t nb_utigs_nucl_mem = 0;
    for (itSeq.first(); !itSeq.isDone(); itSeq.next()) // could be done in parallel, maybe, if we used many unordered_map's with a hash on the query kmer (TODO opt)
    {
        const string& seq = itSeq->toString();
        const string& comment = itSeq->getComment();

        float mean_abundance;
        vector<uint64_t> inc, outc; // incoming and outcoming unitigs
        parse_unitig_header(comment, mean_abundance, inc, outc);

        insert_navigational_vector(incoming,  inc,  utig_counter, incoming_map);
        insert_navigational_vector(outcoming, outc, utig_counter, outcoming_map);

        unitigs.push_back(seq);
        unitigs_mean_abundance.push_back(mean_abundance);

        utig_counter++;
        nb_utigs_nucl += seq.size();
        nb_utigs_nucl_mem += seq.capacity();

        if (seq.size() == kmerSize)
            nb_unitigs_extremities++;
        else
            nb_unitigs_extremities+=2;
    }


    unitigs_traversed.resize(0);
    unitigs_traversed.resize(unitigs.size(), false); // resize "traversed" bitvector, setting it to zero as well

    unitigs_deleted.resize(0);
    unitigs_deleted.resize(unitigs.size(), false); // resize "traversed" bitvector, setting it to zero as well

    // an estimation of memory usage
    uint64_t nb_kmers = unitigs.size();
    uint64_t mem_vec = (unitigs.capacity() * sizeof(string) + nb_utigs_nucl_mem);
    std::cout <<  "Memory usage:" << std::endl;
    std::cout <<  "   " << (sizeof(uint64_t) * incoming.size()) / 1024 / 1024 << " MB keys in incoming dict" << std::endl;
    std::cout <<  "   " << (sizeof(uint64_t) * outcoming.size()) / 1024 / 1024 << " MB keys in outcoming dict" << std::endl;
    std::cout <<  "   " << (sizeof(uint64_t) * incoming_map.size()) / 1024 / 1024 << " MB keys in incoming_map dict" << std::endl;
    std::cout <<  "   " << (sizeof(uint64_t) * outcoming_map.size()) / 1024 / 1024 << " MB keys in outcoming_map dict" << std::endl;
    std::cout <<  "   " <<  mem_vec /1024 /1024 << " MB unitigs nucleotides" << std::endl;
    std::cout <<  "   " <<  (nb_kmers*sizeof(float)) / 1024 / 1024 << " MB unitigs abundances" << std::endl;
    std::cout <<  "   " <<  (2*nb_kmers/8) / 1024 / 1024 << " MB deleted/visited bitvectors" << std::endl;
    std::cout <<  "Estimated total: " <<  (nb_kmers*(sizeof(float) + 2.0/8.0) + sizeof(uint64_t) * ( incoming.size() + outcoming.size() + incoming.size() + outcoming_map.size()) + mem_vec) / 1024 / 1024 << " MB" << std::endl;

    if (nb_utigs_nucl != nb_utigs_nucl_mem)
        std::cout << "unitigs strings size " << nb_utigs_nucl << " vs capacity " << nb_utigs_nucl_mem << std::endl;
}

/*********************************************************************
** METHOD  :
** PURPOSE : creates or completes a graph from parsed command line arguments.
** INPUT   : a bank or a h5 file 
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>::GraphUnitigsTemplate (tools::misc::IProperties* params, bool load_unitigs_after) 
//    : BaseGraph() // call base () constructor // seems to do nothing, maybe it's always called by default
{
   
    string input = params->getStr(STR_URI_INPUT);

    string unitigs_filename = (params->get(STR_URI_OUTPUT) ?
            params->getStr(STR_URI_OUTPUT) :                                                                                                                                                                 System::file().getBaseName (input)
            )+ ".unitigs.fa";          

    // build_visitor_solid has the following defaults:
    // minimizer size of 8. that one is okay
    // the rest needs to be set!
    // accoring to original BCALM2 graph::create string:
    // -in %s -kmer-size %d -minimizer-size %d -mphf none -bloom none -out %s.h5  -abundance-min %d -verbose 1 -minimizer-type %d -repartition-type 1 -max-memory %d %s

    //if ((!params->get(STR_REPARTITION_TYPE)))  // actually this doesn't seem to work, even when repartition-type isn't specified, it's (!params->get()) doesn't return true. so i'm going to force repartition type to 1, as it was in bcalm2
    {
        params->setInt(STR_REPARTITION_TYPE, 1);
        //std::cout << "setting repartition type to 1" << std::endl;;
    }
    //if (!params->get(STR_MINIMIZER_TYPE))
    {
        params->setInt(STR_MINIMIZER_TYPE, 1);
        //std::cout << "setting repartition type to 1" << std::endl;;
    }




    if (system::impl::System::file().getExtension(input) == "h5")
    {
        /* it's not a bank, but rather a h5 file (kmercounted or more), let's complete it to a graph */
        
        if (!System::file().doesExist(input))
            throw system::Exception ("Input file does not exist");

        string output = input;//.substr(0,input.find_last_of(".h5")) + "_new.h5";
        //cout << "To avoid overwriting the input (" << input << "), output will be saved to: "<< output << std::endl;

        cout << "Input is a h5 file (we assume that it contains at least the solid kmers).\n"; 
        
        /** We create a storage instance. */
        /* (this is actually loading, not creating, the storage at "uri") */
        BaseGraph::setStorage (StorageFactory(BaseGraph::_storageMode).create (output , false, false));
    
        /** We get some properties. */
        BaseGraph::_state     = (typename GraphUnitigsTemplate<span>::StateMask) atol (BaseGraph::getGroup().getProperty ("state").c_str());
        
        BaseGraph::_kmerSize  = atol (BaseGraph::getGroup().getProperty ("kmer_size").c_str());

        if (BaseGraph::_kmerSize == 0) /* try the dsk group -> maybe it's a dsk h5 file, not a minia one */
            BaseGraph::_kmerSize  =    atol (BaseGraph::getGroup("dsk").getProperty ("kmer_size").c_str());
        
        modelK = new Model(BaseGraph::_kmerSize);
        modelKdirect= new ModelDirect(BaseGraph::_kmerSize);

        // also assume kmer counting is done
        setState(GraphUnitigsTemplate<span>::STATE_SORTING_COUNT_DONE);
        
        /** We get library information in the root of the storage. */
        string xmlString = BaseGraph::getGroup().getProperty ("xml");
        stringstream ss; ss << xmlString;   IProperties* props = new Properties(); LOCAL(props);
        props->readXML (ss);  BaseGraph::getInfo().add (1, props);
        
        /** We configure the data variant according to the provided kmer size. */
        BaseGraph::setVariant (BaseGraph::_variant, BaseGraph::_kmerSize);

        /* call the configure visitor to load everything (e.g. solid kmers, MPHF, etc..) that's been done so far */
        boost::apply_visitor ( configure_visitor<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>>(*this, BaseGraph::getStorage()),  *(GraphDataVariantFast<span>*)BaseGraph::_variant);

        build_unitigs_postsolid(unitigs_filename, params);
      
        if (load_unitigs_after) 
            load_unitigs(unitigs_filename);


    }
    else
    {
        /** We get the kmer size from the user parameters. */
        BaseGraph::_kmerSize = params->getInt (STR_KMER_SIZE);
        modelK = new Model(BaseGraph::_kmerSize);
        modelKdirect= new ModelDirect(BaseGraph::_kmerSize);
        size_t integerPrecision = params->getInt (STR_INTEGER_PRECISION);

        /** We configure the data variant according to the provided kmer size. */
        BaseGraph::setVariant (BaseGraph::_variant, BaseGraph::_kmerSize, integerPrecision);

        /** We build a Bank instance for the provided reads uri. */
        bank::IBank* bank = Bank::open (params->getStr(STR_URI_INPUT));

        /** We build the graph according to the wanted precision. */
        boost::apply_visitor ( build_visitor_solid<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>>(*this, bank,params),  *(GraphDataVariantFast<span>*)BaseGraph::_variant);

        build_unitigs_postsolid(unitigs_filename, params);

        if (load_unitigs_after)
            load_unitigs(unitigs_filename);
    }
}

template<size_t span>
GraphUnitigsTemplate<span>::GraphUnitigsTemplate ()
    : GraphTemplate<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>>()
{
}

template<size_t span>
GraphUnitigsTemplate<span>::GraphUnitigsTemplate (const GraphUnitigsTemplate<span>& graph)
    : GraphTemplate<NodeFast<span>,EdgeFast<span>,GraphDataVariantFast<span>>(graph)
{
    // will call Graph's constructor
    std::cout << "GraphU copy-constructor called" << std::endl;

}

/*********************************************************************
** PURPOSE : copy assignment operator
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>& GraphUnitigsTemplate<span>::operator= (GraphUnitigsTemplate<span> const& graph)
{
    std::cout <<"assignment constructor called" << std::endl;
    if (this != &graph)
    {
        BaseGraph::_kmerSize        = graph._kmerSize;
        BaseGraph::_storageMode     = graph._storageMode;
        BaseGraph::_name            = graph._name;
        BaseGraph::_info            = graph._info;
        BaseGraph::_mphfKind        = graph._mphfKind;
        BaseGraph::_state           = graph._state;

        BaseGraph::setStorage (graph._storage);

        if (graph._variant)  {  *((GraphDataVariantFast<span>*)BaseGraph::_variant) = *((GraphDataVariantFast<span>*)graph._variant);  }

        // don't forget those!
        // I garantee that bugs will occur if i add a GraphUnitigs member variable and forget to copy it here
        // DO ALSO THE MOVE FUNCTION BELOW
    
        incoming  = graph.incoming;
        outcoming = graph.outcoming;
        incoming_map  = graph.incoming_map;
        outcoming_map = graph.outcoming_map;
        unitigs = graph.unitigs;
        unitigs_mean_abundance = graph.unitigs_mean_abundance;
        unitigs_traversed = graph.unitigs_traversed;
        unitigs_deleted = graph.unitigs_deleted;
        modelK = graph.modelK;
        modelKdirect = graph.modelKdirect;
        nb_unitigs = graph.nb_unitigs;
        nb_unitigs_extremities = graph.nb_unitigs_extremities;
        
    }
    return *this;
}

/*********************************************************************
** PURPOSE : move assignment operator
** this seems important, if it was not there, in Minia, the line "graph = create (..)" would incur a copy of the graph :/ didn't investigate further, implementing that move operator was enough to prevent the copy.
// also the copy-and-swap idiom didn't seem to help, i got infinite loop because of effect described here http://stackoverflow.com/questions/25942131/should-copy-assignment-operator-leverage-stdswap-as-a-general-rule
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>& GraphUnitigsTemplate<span>::operator= (GraphUnitigsTemplate<span> && graph)
{
    std::cout <<"move constructor called" << std::endl;
    if (this != &graph)
    {
        BaseGraph::_kmerSize        = graph._kmerSize;
        BaseGraph::_storageMode     = graph._storageMode;
        BaseGraph::_name            = graph._name;
        BaseGraph::_info            = graph._info;
        BaseGraph::_mphfKind        = graph._mphfKind;
        BaseGraph::_state           = graph._state;

        BaseGraph::setStorage (graph._storage);

        if (graph._variant)  {  *((GraphDataVariantFast<span>*)BaseGraph::_variant) = *((GraphDataVariantFast<span>*)graph._variant);  }

        // don't forget those!
        // I garantee that bugs will occur if i add a GraphUnitigs member variable and forget to copy it here
    
        incoming  = std::move(graph.incoming);
        outcoming = std::move(graph.outcoming);
        incoming_map  = std::move(graph.incoming_map);
        outcoming_map = std::move(graph.outcoming_map);
        unitigs = std::move(graph.unitigs);
        unitigs_mean_abundance = std::move(graph.unitigs_mean_abundance);
        unitigs_traversed = std::move(graph.unitigs_traversed);
        unitigs_deleted = std::move(graph.unitigs_deleted);
        modelK = std::move(graph.modelK);
        modelKdirect = std::move(graph.modelKdirect);
        nb_unitigs = std::move(graph.nb_unitigs);
        nb_unitigs_extremities = std::move(graph.nb_unitigs_extremities);
        
    }
    return *this;
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
GraphUnitigsTemplate<span>::~GraphUnitigsTemplate<span> ()
{
    // base deleter already called
    std::cout <<"unitigs graph destructor called" << std::endl;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphUnitigsTemplate<span>::remove ()
{
    std::cout << "GraphU remove called" << std::endl;
    BaseGraph::getStorage().remove();
}


    /* just a note: nothing weird with having 1-in, 1-out nodes in GraphUnitigs: consider that example:
     *        --------
     *                 \
     *                  v
     *  -------[node] -> -----------
     *
     * [node] is 1-in and 1-out yet compactions were fine, it's just that there is in-branching in the following node.
     */

template<size_t span>
GraphVector<EdgeGU> GraphUnitigsTemplate<span>::getEdges (const NodeGU& source, Direction direction)  const
{
    bool debug = false;

    if (debug)
    {
        std::cout << "graphU getEdges called, on source: " << toString(source) << " unitig: " << source.unitig << " pos: " << (source.pos==UNITIG_BEGIN?"beg":"end") << " strand: " << source.strand << " dir " << direction << std::endl;
    }

    if (source.pos == UNITIG_INSIDE)
    {
        std::cout << "Error: called getEdges on a node that's UNITIG_INSIDE: " << toString(source) << std::endl; exit(1);
    }
 
    GraphVector<EdgeGU> res;

    res.resize(0);
    
    unsigned int kmerSize = BaseGraph::_kmerSize;
    unsigned int seqSize = unitigs[source.unitig].size();
    
    bool same_orientation = node_in_same_orientation_as_in_unitig(source);
    bool pos_begin = source.pos & UNITIG_BEGIN;
    bool pos_end = source.pos & UNITIG_END;

    if ((!pos_begin) && (!pos_end))
    {
        std::cout << "weird node position: " << source.pos << " unitig length: " << seqSize;
        exit(1);
    }

    if ((unsigned int)seqSize > (unsigned int)kmerSize)
    {
        // unitig: [kmer]-------
        if (same_orientation && (direction & DIR_OUTCOMING) && pos_begin) 
        {
            Unitig_pos pos = UNITIG_INSIDE;
            if (seqSize == kmerSize + 1) pos = UNITIG_END;
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.unitig, source.pos, source.strand, source.unitig, pos, source.strand, DIR_OUTCOMING);
            if (debug) std::cout << "found success of [kmer]---" << std::endl;
        }

        // unitig: [kmer rc]-------
        if ((!same_orientation) && (direction & DIR_INCOMING) && pos_begin) 
        {
            Unitig_pos pos = UNITIG_INSIDE;
            if (seqSize == kmerSize + 1) pos = UNITIG_END;
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.unitig, source.pos, source.strand, source.unitig, pos, source.strand, DIR_INCOMING); // not sure about the dest.strand in those cases, so i'm setting to source.strand, we'll see. (applies to all 3 other cases) It doesn't matter in Minia anyway, we don't use nodes inside unitigs
            if (debug) std::cout << "found success of [kmer rc]---" << std::endl;
        }

        // unitig: ----------[kmer]
        if ((same_orientation) && (direction & DIR_INCOMING) && pos_end)
        {
            Unitig_pos pos = UNITIG_INSIDE;
            if (seqSize == kmerSize + 1) pos = UNITIG_BEGIN;
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.unitig, source.pos, source.strand, source.unitig, pos, source.strand, DIR_INCOMING);
            if (debug) std::cout << "found predec of --------[kmer]" << std::endl;
        }

        // unitig: ----------[kmer rc]
        if ((!same_orientation) && (direction & DIR_OUTCOMING) && pos_end)
        {
            Unitig_pos pos = UNITIG_INSIDE;
            if (seqSize == kmerSize + 1) pos = UNITIG_BEGIN;
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.unitig, source.pos, source.strand, source.unitig, pos, source.strand, DIR_OUTCOMING);
            if (debug) std::cout << "found predec of --------[kmer rc]" << std::endl;
        }
    }
    else
    {
        pos_begin = pos_end = true; // necessary fix that sohuld have happened in parse_unitig_header
    }
    
    // otherwise, that extremity kmer has neighbors at are also extremities.
    // so, mutate to get all 4 outneighrs, and test for their existence in the utigs_map
    
    auto functor = [&](range<std::vector<uint32_t>::const_iterator >&& edges, Direction dir)
    {
        auto it = edges.begin();
        if (it == edges.end()) return;
        for (; it != edges.end(); it++)
        {
            auto edge_packed = *it;
            ExtremityInfo li(edge_packed);

            if (unitigs_deleted[li.unitig]) 
            {
                if (debug)
                    std::cout << "found deleted neighbor unitig "<<  li.unitig <<" (kmer: " << this->/*not putting this this crashes gcc 4.7 */toString(NodeGU(li.unitig, li.pos)) << ")" << std::endl;
                continue;
            }

            uint64_t unitig = li.unitig;
            Unitig_pos pos = li.pos;
            bool rc = li.rc;

            if (!same_orientation)
                rc = !rc;

            kmer::Strand strand = rc?STRAND_REVCOMP:STRAND_FORWARD;

            if (debug) 
            {
                NodeGU node(unitig, pos, strand);
                std::cout << "[out-of-unitig getEdges], found neighbor " << this->/*not putting this this crashes gcc 4.7 */toString(node) << " dir " << dir << std::endl;
            }
    
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.unitig, source.pos, source.strand, unitig, pos, strand, dir);
        }
    }; 

    // TODO: write an unit test for the bug where not putting paranthesis there: pos_end && (...) was still working, but causes more edges to be returned than necessary. (only in DIR_END i think)
    if (pos_end && (((direction & DIR_OUTCOMING) && same_orientation) || ( (direction & DIR_INCOMING) && (!same_orientation) ) ))
    {
        // nodes to the right of a unitig (outcoming)
        Direction dir = same_orientation?DIR_OUTCOMING:DIR_INCOMING;
        functor(get_from_navigational_vector(outcoming, source.unitig, outcoming_map), dir);
    }
    if (pos_begin && (((direction & DIR_INCOMING) && same_orientation) || ( (!same_orientation) && (direction & DIR_OUTCOMING)) ))
    {
        // nodes to the left of a unitig (incoming)
        Direction dir = same_orientation?DIR_INCOMING:DIR_OUTCOMING;
        functor(get_from_navigational_vector(incoming, source.unitig, incoming_map), dir);
    }

    // sanity check on output, due to limitation on GraphVector nmber of elements
    // TODO address that. i don't like the potential performance hit here.
    if (res.size() > 16)
    { 
        std::cout << "Error : more than 16 edges (" << res.size() << ") out of node, not supported (already more than 8 is strange)" << std::endl; 
        std::cout << "graphU getEdges was called on source: " << toString(source) << " unitig: " << source.unitig << " pos: " << (source.pos==UNITIG_BEGIN?"beg":"end") << " strand: " << source.strand << " dir " << direction << std::endl;
        exit(1);
    }
    
    return res;
}

/* this function isn't the most efficient, but then again, Minia doesn't use it */
// I think I only implemented it so that it passes simple tests (TestDebruijnUnitigs.cpp)
template<size_t span>
GraphVector<NodeGU> GraphUnitigsTemplate<span>::getNodes (const NodeGU &source, Direction direction)  const
{
    GraphVector<NodeGU> nodes;
    GraphVector<EdgeGU> edges = getEdges (source, direction);
    nodes.resize(edges.size());
    for (unsigned int i = 0; i < edges.size(); i++)
    {
        nodes[i] = edges[i].to;
    }
    return nodes;
}


template<size_t span>
unsigned char GraphUnitigsTemplate<span>::countNeighbors (const NodeGU &source, Direction direction)  const
{
    // for the sake of no duplication and removing bugs, i'm de-optimizing this function for now.
    GraphVector<EdgeGU> edges = getEdges(source, direction);
    return edges.size();
}

template<size_t span>
void GraphUnitigsTemplate<span>::countNeighbors (const NodeGU &source, size_t &in, size_t &out)  const
{
    std::cout << "GraphU countNeighbors source,in,out not implememented" << std::endl;exit(1);
}

template<size_t span>
NodeGU GraphUnitigsTemplate<span>::getNode (const NodeGU& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{
    std::cout << "GraphU getNode source,dir,nt,exists  not implememented" << std::endl;exit(1);
    return NodeGU();
}

template<size_t span>
GraphIterator<NodeGU> GraphUnitigsTemplate<span>::getNodes () const
{
    /* emulates iteration of nodes à la original GATB Graph */
    /* except that here, we only iterate the extremities of unitigs */
    class NodeIterator : public tools::dp::ISmartIterator<NodeGU>
    {
        public:
            NodeIterator (const std::vector<std::string>& unitigs /* just to get lengths*/, const std::vector<bool>& unitigs_deleted, unsigned int k, unsigned int nb_unitigs_extremities) 
                :  _nbItems(nb_unitigs_extremities), _rank(0), _isDone(true), unitigs(unitigs), unitigs_deleted(unitigs_deleted), k(k), nb_unitigs(unitigs.size()) {  
                    this->_item->strand = STRAND_FORWARD;  // iterated nodes are always in forward strand.
                }

            ~NodeIterator ()  {  }

            u_int64_t rank () const { return _rank; }

            void update_item()
            {
                this->_rank ++;
                this->_item->unitig = it/2;
                this->_item->pos = (it&1)?UNITIG_END:UNITIG_BEGIN;
            }

            /** \copydoc  Iterator::first */
            void first()
            {
                it = 0;
                while (unitigs_deleted[it/2] && it < 2*nb_unitigs) it++;
                _rank   = 0;
                _isDone = it >= (2*nb_unitigs);

                if (!_isDone)
                    update_item();
            }

            /** \copydoc  Iterator::next */
            void next()
            {
                do
                {
                    it++;
                    if ((it < 2*nb_unitigs) && unitigs[it/2].size() == k) // takes care of the case where the unitig is just a kmer
                        it++;
                } while ((it < 2*nb_unitigs) && unitigs_deleted[it/2]);
                _isDone = it >= (2*nb_unitigs);
                if (!_isDone)
                    update_item();
            }

            /** \copydoc  Iterator::isDone */
            bool isDone() { return _isDone;  }

            /** \copydoc  Iterator::item */
            NodeGU& item ()  {  return *(this->_item);  }

            void setItem (NodeGU& i)
            {
                /** We set the node item to be set for the current iterator. */
                this->_item = &i;
                this->_item->strand = STRAND_FORWARD;
            }

            /** */
            u_int64_t size () const { return _nbItems; }

        private:
            uint64_t it;
            u_int64_t _nbItems;
            u_int64_t _rank;
            bool      _isDone;
            const std::vector<std::string>& unitigs;
            const std::vector<bool>& unitigs_deleted;
            unsigned int k;
            unsigned int nb_unitigs;
    };

    return new NodeIterator(unitigs, unitigs_deleted, BaseGraph::_kmerSize, nb_unitigs_extremities);
}

template<size_t span> 
bool GraphUnitigsTemplate<span>::isNodeDeleted(const NodeGU& node) const
{
    return unitigs_deleted[node.unitig];
}



// emulation of MPHF, but it's not used in minia
template<size_t span> 
unsigned long GraphUnitigsTemplate<span>::nodeMPHFIndex(const NodeGU& node) const 
{
    return (((node.pos == UNITIG_BEGIN )) ? 0 : 1) + (node.unitig << 1); 
}

/* warning: will delete the whole simple path */
template<size_t span>
void GraphUnitigsTemplate<span>::deleteNode (NodeGU& node) 
{
    unitigDelete (node);
}

template<size_t span>
void GraphUnitigsTemplate<span>::cacheNonSimpleNodes(unsigned int nbCores, bool verbose) 
{
    // nothing to do with unitigs flavor.
}
 
template<size_t span>
void GraphUnitigsTemplate<span>::deleteNodesByIndex(vector<bool> &bitmap, int nbCores, gatb::core::system::ISynchronizer* synchro) const
{
    std::cout << "deleteNodesByIndex called, shouldn't be." << std::endl; 
    exit(1);
}

/********************************************************************************/
template<size_t span>
bool GraphUnitigsTemplate<span>::contains (const NodeGU& item) const
{
    std::cout << "contains() not implemeneted in GraphUnitigs" << std::endl; exit(1);
    return false;
}

template<size_t span>
bool GraphUnitigsTemplate<span>::
node_in_same_orientation_as_in_unitig(const NodeGU& node) const
{
    return node.strand == kmer::STRAND_FORWARD;
}

template<size_t span>
std::string GraphUnitigsTemplate<span>::toString (const NodeGU& node) const
{
    const std::string& seq = unitigs[node.unitig];
    int kmerSize = BaseGraph::_kmerSize;

    if (node.pos == UNITIG_INSIDE)
    {    return "[GraphUnitigs.toString cannot print an UNITIG_INSIDE]"; }
    
    string node_str;
    if (node.pos & UNITIG_BEGIN)
        node_str = seq.substr(0,kmerSize);
    else
        node_str = seq.substr(seq.size()-kmerSize);

    if (node.strand != kmer::STRAND_FORWARD)
        node_str = revcomp(node_str);

    return node_str;
}

// high-level functions that used to be in Simplifications.cpp
template<size_t span>
bool GraphUnitigsTemplate<span>::
isLastNode(const NodeGU& node, Direction dir) const
{
    if (unitigs[node.unitig].size() == BaseGraph::_kmerSize) // special case.
        return true;

    bool same_orientation = node_in_same_orientation_as_in_unitig(node);
    Unitig_pos pos = node.pos;

    // cases where, following that unitig direction, we're already at the last node
    if ((same_orientation    && (pos & UNITIG_END) && dir == DIR_OUTCOMING) ||
        (same_orientation    && (pos & UNITIG_BEGIN) && dir == DIR_INCOMING) ||
        ((!same_orientation) && (pos & UNITIG_END) && dir == DIR_INCOMING) ||
        ((!same_orientation) && (pos & UNITIG_BEGIN) && dir == DIR_OUTCOMING))
        return true;
    return false;
}

template<size_t span>
bool GraphUnitigsTemplate<span>::
isFirstNode(const NodeGU& node, Direction dir) const
{
    // special case
    if (unitigs[node.unitig].size() == BaseGraph::_kmerSize)
    {
        return true;
    }

    return !isLastNode(node,dir);
}

template<size_t span>
double GraphUnitigsTemplate<span>::
simplePathMeanAbundance     (const NodeGU& node, Direction dir) 
{
    if (isLastNode(node,dir))
    {
        if (!isFirstNode(node,dir))
            return 0;
        else // single-k-mer unitig
            return unitigMeanAbundance(node);
        }

    float coverage = 0;
    int endDegree;
    int seqLength = 0;
    simplePathLongest_avance(node, dir, seqLength, endDegree, false /*markDuringTraversal*/, coverage);
    return coverage / (float)seqLength;
}

template<size_t span>
double GraphUnitigsTemplate<span>::
unitigMeanAbundance     (const NodeGU& node) const
{
    return unitigs_mean_abundance[node.unitig];
}

template<size_t span>
unsigned int GraphUnitigsTemplate<span>::
simplePathLength            (const NodeGU& node, Direction dir) 
{
    float coverage = 0;
    int endDegree;
    int seqLength = 0;
    simplePathLongest_avance(node, dir, seqLength, endDegree, false /*markDuringTraversal*/, coverage);
    return seqLength;
}

template<size_t span>
unsigned int GraphUnitigsTemplate<span>::
unitigLength            (const NodeGU& node, Direction dir) const
{
    const std::string& seq = unitigs[node.unitig];
    bool same_orientation = node_in_same_orientation_as_in_unitig(node);

    int length;
    if (isLastNode(node,dir))
        length = 0;
    else
        length = unitigs[node.unitig].size() - BaseGraph::_kmerSize;
    return length;
}

template<size_t span>
NodeGU GraphUnitigsTemplate<span>::
unitigLastNode          (const NodeGU& node, Direction dir) const
{
    //const std::string& seq = unitigs[node.unitig];

    //std::cout << "lastnode" << toString(node) << " dir " << dir  << std::endl;
    
    if (isLastNode(node,dir))
        return node;

    // otherwise we know that we just take the other extremity, easy peasy
    
    NodeGU res;
    if ((node.pos & UNITIG_BEGIN))
        res = NodeGU(node.unitig, UNITIG_END, node.strand);
    else
        res = NodeGU(node.unitig, UNITIG_BEGIN, node.strand);
    return res;
}

template<size_t span>
NodeGU GraphUnitigsTemplate<span>::
simplePathLastNode          (const NodeGU& node, Direction dir) 
{
    std::vector<NodeGU> nodesList;
    int seqLength = 0, endDegree;
    float coverage = 0;
    simplePathLongest_avance(node, dir, seqLength, endDegree, false /*markDuringTraversal*/, coverage, nullptr, &nodesList);
    if (nodesList.size() == 0)
    {
        assert(isLastNode(node,dir));
        return node;
    }
    return nodesList.back();
}


template<size_t span>
void GraphUnitigsTemplate<span>::
unitigDelete (NodeGU& node) 
{
    unitigs_deleted[node.unitig] = true;
    //std::cout << "GraphU deleted unitig " << node.unitig << " seq: "  << unitigs[node.unitig] << std::endl; 
}


// wrapper that only puts the node into nodesDeleter
template<size_t span>
void GraphUnitigsTemplate<span>::
unitigDelete (NodeGU& node, Direction dir, NodesDeleter<NodeGU, EdgeGU, GraphUnitigsTemplate<span>>& nodesDeleter) 
{
    nodesDeleter.onlyListMethod = true; // a bit inefficient to always tell the deleter to be in that mode, but so be it for now. just 1 instruction, won't hurt.
    //std::cout << "GraphU queuing to delete unitig " << node.unitig << " seq: "  << unitigs[node.unitig] << " mean abundance " << unitigMeanAbundance(node) << std::endl; 
    nodesDeleter.markToDelete(node);
}

template<size_t span>
void GraphUnitigsTemplate<span>::
simplePathDelete (NodeGU& node, Direction dir, NodesDeleter<NodeGU, EdgeGU, GraphUnitigsTemplate<span>>& nodesDeleter) 
{
    std::vector<NodeGU> nodesList;
    int seqLength = 0, endDegree;
    float coverage = 0;
    unitigDelete(node, dir, nodesDeleter);
    simplePathLongest_avance(node, dir, seqLength, endDegree, false /*markDuringTraversal*/, coverage, nullptr, &nodesList);
    for (auto cur_node : nodesList)
    {
        unitigDelete(cur_node, dir, nodesDeleter);
    }
}

/* actually I don't think this function is called at all */
template<size_t span>
std::string GraphUnitigsTemplate<span>::
unitigSequence (const NodeGU& node, bool& isolatedLeft, bool& isolatedRight) const
{
    const string& seq = unitigs[node.unitig];

    //std::cout << " seq " << seq << " node " << toString(node) << std::endl;
    int kmerSize = BaseGraph::_kmerSize;
    NodeGU left = NodeGU(node.unitig, UNITIG_BEGIN);
    NodeGU right = NodeGU(node.unitig, UNITIG_END);

    isolatedLeft  = (indegree(left)   == 0);
    isolatedRight = (outdegree(right) == 0);
    return seq;
}

// keep traversing unitigs as far as we can.
/* arguments:
 * starting node 
 * output sequence length (will be reset if previously set)
 * output end degree
 * whether to mark during traversal
 * total coverage of simple path (NOT normalized!)
 * output sequence pointer (could be NULL if we don't care about it)
 * all nodes extremities of simple paths traversed
 */
template<size_t span>
void GraphUnitigsTemplate<span>::
simplePathLongest_avance(const NodeGU& node, Direction dir, int& seqLength, int& endDegree, bool markDuringTraversal, float& coverage, string* seq, std::vector<NodeGU> *nodesList) 
{
    bool debug = false;
    if (debug)
        std::cout << "simplePathLongest_avance called, node " << toString(node) << " dir " << dir << std::endl;

    seqLength = 0;
    unsigned int kmerSize = BaseGraph::_kmerSize;
    NodeGU cur_node = node;
    if (!isLastNode(cur_node,dir))
    {
        // first node in unitig may have in-branching, it's fine for a simple path traversal. we'll just go to last node and record the sequence of that unitig
        int unitigLength = unitigs[node.unitig].size();

        bool same_orientation = node_in_same_orientation_as_in_unitig(node);

        if (seq != nullptr)
        {
            string new_seq = unitigs[node.unitig];

            if (!same_orientation)
                new_seq = revcomp(new_seq);

            if (dir == DIR_OUTCOMING)
                *seq += new_seq.substr(kmerSize-1);
            else
                *seq = new_seq.substr(0,new_seq.size()-(kmerSize-1)) + *seq;
        }

        // length is all the kmers of that unitig, except first one
        seqLength += unitigLength - kmerSize;
        // for coverage: it includes the abundance of the first kmer, because there's no way to exclude it, hence the +1
        coverage += unitigMeanAbundance(cur_node) * (unitigLength - kmerSize + 1);
       
        if (debug)
            std::cout << "simplePathLongest_avance was at a first node = " << toString(cur_node) << " strand " << cur_node.strand << " so traversed unitig of length " << unitigLength << std::endl;

        cur_node = unitigLastNode(node,dir);       
        
        if (debug)
            std::cout << "simplePathLongest_avance now at last node = " << toString(cur_node) << " strand " << cur_node.strand << std::endl;

        if (nodesList != nullptr)
            nodesList->push_back(cur_node);
    }

    if (markDuringTraversal) // no need to mark, that unitig should already be marked by minia, but doing it anyway just to be safe
        unitigMark(cur_node);

    set<uint64_t> traversed_unitigs;
    
    while (true)
    { // invariant here: cur_node is the last node of a unitig

        assert(isLastNode(cur_node,dir));

        GraphVector<EdgeGU> neighbors = this->neighborsEdge (cur_node, dir);
        endDegree = neighbors.size();
        
        /** We check we have no outbranching. */
        if (endDegree != 1)
        {
            if (debug)
                std:: cout << "simplePathLongest_avance stopped because # neighbor of node " << toString(cur_node) << ": " << neighbors.size() << std::endl;
            return;
        }
      
        uint64_t neighbor_unitig = neighbors[0].to.unitig;
        if (traversed_unitigs.find(neighbor_unitig) != traversed_unitigs.end())
        {
            if (debug)
                std::cout << "simplePathLongest_avance loop" << std::endl;
            break;
        }
        traversed_unitigs.insert(neighbor_unitig);

        bool same_orientation = node_in_same_orientation_as_in_unitig(neighbors[0].to);
            
        int unitigLength = unitigs[neighbor_unitig].size();

        if (debug)
            std::cout << "simplePathLongest_avance continues now at a last node = " << toString(cur_node) << " strand " << cur_node.strand << " of unitig " << cur_node.unitig << " length " << unitigs[cur_node.unitig].size() << ", neighbor.to " << toString(neighbors[0].to) << " strand " << neighbors[0].to.strand << " of unitig " << neighbors[0].to.unitig << " new seq length: " << unitigLength << std::endl;

        // fix for 1-bit encoded unitig position. That fix could have happened in unitig_parse_header but i didn't want to encode pos in 2 bits. Also, could have happened in NodeGU constructor, but didn't want to waste time checking unitig size there
        if (unitigs[neighbors[0].to.unitig].size() == kmerSize)
            neighbors[0].to.pos = UNITIG_BOTH;

        int npos = neighbors[0].to.pos;

        // some sanity checks      
        if (neighbors[0].to.pos != UNITIG_BOTH) 
        {
            // FIXME: there is a bcalm bug. see strange_seq.fa. i'll send it to antoine, but meanwhile, i'm coding a workaround
            if (0)
            {
            if (dir==DIR_INCOMING )
                assert( (npos == UNITIG_END   && same_orientation) || (npos == UNITIG_BEGIN && (!same_orientation)));
            else
                assert( (npos == UNITIG_BEGIN && same_orientation) || (npos == UNITIG_END   && (!same_orientation)));
            }
            if (dir==DIR_INCOMING )
            {
                if (!( (npos == UNITIG_END   && same_orientation) || (npos == UNITIG_BEGIN && (!same_orientation))))
                {
                    unitigDelete(neighbors[0].to);
                    return;
                }
            }
            else
            {
                if (!( (npos == UNITIG_BEGIN && same_orientation) || (npos == UNITIG_END   && (!same_orientation))))
                {
                    unitigDelete(neighbors[0].to);
                    return;
                }
            }

        }
       
        GraphVector<EdgeGU> in_neighbors_vec = this->neighborsEdge (neighbors[0].to, reverse(dir));
        int in_neighbors = in_neighbors_vec.size();

        assert(in_neighbors >= 1);
        /** We check we have no in-branching. */
        if (in_neighbors > 1) 
        {
            if (debug)
                std:: cout << "simplePathLongest_avance stopped at " << toString(cur_node) << " because of in-branching " << in_neighbors << std::endl;
            return;
        } 

        NodeGU last_node = unitigLastNode(neighbors[0].to, dir);
        cur_node = last_node;

        if (nodesList != nullptr)
            nodesList->push_back(cur_node);

        // append the sequence (except the overlap part, of length k-1.
        if (seq != nullptr)
        {
            string new_seq = unitigs[cur_node.unitig];
            if (!same_orientation)
                new_seq = revcomp(new_seq);

            if (dir == DIR_OUTCOMING)
                *seq += new_seq.substr(kmerSize-1);
            else
                *seq = new_seq.substr(0,new_seq.size()-(kmerSize-1)) + *seq;
        }

        seqLength += unitigLength - (kmerSize-1);
        coverage += unitigMeanAbundance(cur_node) * (unitigLength - kmerSize + 1); // here too, coverage is computed according to whole unitig

        if (markDuringTraversal&& unitigIsMarked(cur_node)) // just a debug, can be removed
        {
            //std::cout << "marked node during a simple path traversal, that shouldn't happen. Maybe it's a perfect loop." << std::endl;
            return;
        }

        if (markDuringTraversal)
            unitigMark(cur_node);
    }
}

/* returns the longest simple path; may have to traverse multiple unitigs, due to some branches being deleted */
template<size_t span>
std::string GraphUnitigsTemplate<span>::
simplePathBothDirections(const NodeGU& node, bool& isolatedLeft, bool& isolatedRight, bool markDuringTraversal, float &coverage) 
{
    string seq = unitigs[node.unitig];
    
    int kmerSize = BaseGraph::_kmerSize;
    float midTotalCoverage = unitigMeanAbundance(node) * (seq.size() - kmerSize + 1);

    NodeGU left(node.unitig, UNITIG_BEGIN);
    NodeGU right(node.unitig, UNITIG_END);

    //std::cout << "starting seq " << seq << " (from node " << toString(node) << ") left: " << toString(left) << " right: " << toString(right) << std::endl;

    if (markDuringTraversal)
        unitigMark(left);

    string seqRight = "", seqLeft = "";
    int endDegreeLeft, endDegreeRight;
    float rightTotalCoverage = 0, leftTotalCoverage = 0;
    int lenSeqRight = 0, lenSeqLeft = 0;
    simplePathLongest_avance (right, DIR_OUTCOMING, lenSeqRight, endDegreeRight, markDuringTraversal, rightTotalCoverage, &seqRight);
    simplePathLongest_avance (left, DIR_INCOMING, lenSeqLeft, endDegreeLeft, markDuringTraversal, leftTotalCoverage, &seqLeft);

    isolatedLeft  = (endDegreeLeft == 0);
    isolatedRight = (endDegreeRight == 0);

    // glue everything together
    seq = seqLeft + seq + seqRight;
    coverage = (rightTotalCoverage + leftTotalCoverage + midTotalCoverage) / (seq.size() - kmerSize + 1);
    return seq;
}


// used to flag simple path as traversed, in minia
// marks the whole unitig, not just an extremity
template<size_t span>
void GraphUnitigsTemplate<span>::
unitigMark            (const NodeGU& node) 
{
    unitigs_traversed[node.unitig] = true;
} 

template<size_t span>
bool GraphUnitigsTemplate<span>::
unitigIsMarked        (const NodeGU& node) const 
{
    return unitigs_traversed[node.unitig];
}

/* when debugging goes wrong;. print the whole graph. */
template<size_t span>
void GraphUnitigsTemplate<span>::
debugPrintAllUnitigs() const
{
    std::cout << "Debug: printing all graph unitigs and status" << std::endl;
    for (unsigned int i = 0; i < nb_unitigs; i++)
    {
        std::cout << "unitig " << i << " (length: " << unitigs[i].size() << ") " << (unitigs_deleted[i]?"[deleted]":"") << " links: ";


        for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
        { 
            GraphVector<EdgeGU> neighbors;
            if (dir == DIR_OUTCOMING)
            {
                NodeGU cur_node(i, UNITIG_END);
                neighbors = this->neighborsEdge (cur_node, dir);
                std::cout << "out: ";
            }
            else
            {
                NodeGU cur_node(i, UNITIG_BEGIN);
                neighbors = this->neighborsEdge (cur_node, dir);
                std::cout << "in: ";
            }
            for (size_t i = 0; i < neighbors.size(); i++)
            {
                string pos = (neighbors[i].to.pos == UNITIG_BEGIN)?"beg":((neighbors[i].to.pos == UNITIG_END)?"end":"other");
                std::cout << neighbors[i].to.unitig << "(" << pos << ") ";
            }
            std::cout <<  "    ";
        }
        std::cout << std::endl;
    }
    std::cout << "Done printing graph unitigs" << std::endl;
}

/* slow, O(|unitigs|) construction of a node given its kmer. for debug only */
template<size_t span>
NodeGU GraphUnitigsTemplate<span>::
debugBuildNode(string startKmer) const
{
    for (unsigned int i = 0; i < unitigs.size(); i++)
    {
        string unitig = unitigs[i];
        for (int rc = 0; rc < 2; rc++)
        {
            if (rc == 1) 
                unitig = revcomp(unitig);
            if (unitig.substr(0, BaseGraph::_kmerSize) == startKmer)
            {
                return NodeGU(i,rc?UNITIG_END:UNITIG_BEGIN, rc?STRAND_REVCOMP:STRAND_FORWARD);
            }
            if (unitig.substr(unitig.size() - BaseGraph::_kmerSize) == startKmer)
            {
                return NodeGU(i,rc?UNITIG_BEGIN:UNITIG_END, rc?STRAND_REVCOMP:STRAND_FORWARD);
            }
        }
    }
    std::cout << "could not build node correspoding to kmer" << startKmer << ", it's not an unitig end." << std::endl;
    exit(1);
}


/*
 *
 *
 *
 *
 *
 *
 * boilerplate code that's common with Graph.hpp but i'm duplicating because too lazy to do a proper abstract class for GraphAbstract (and would be tedious no)
 in other words: would be good to factorize with Graph.hpp but i think it'd require having a GraphAbstract and I'm not ready for that kind of design pattern yet.
 *
 * no algorithmic changes needed in GraphUnitigs
 *
 */


template<size_t span>
bool GraphUnitigsTemplate<span>::isBranching (const NodeGU& node) const
/** REMARKS : it's a bit of a misnomer, but in GATB language, Branching means "anything not simple". I.e. simple means outdegree=indegree=1 and branching is every other degree combination, including 0's */
{
    size_t in, out;
    degree(node, in, out); 
    return (! (in==1 && out==1));
}

template<size_t span>
bool GraphUnitigsTemplate<span>::isSimple (EdgeGU& edge) const
{    return this->outdegree(edge.from)==1  &&  this->indegree(edge.to)==1;}

template<size_t span>
size_t GraphUnitigsTemplate<span>::indegree (const NodeGU& node) const  {  return degree(node, DIR_INCOMING);   }

template<size_t span>
size_t GraphUnitigsTemplate<span>::outdegree (const NodeGU& node) const  {  return degree(node, DIR_OUTCOMING);  }

template<size_t span>
size_t GraphUnitigsTemplate<span>::degree (const NodeGU& node, Direction dir) const  {  return countNeighbors(node, dir);  } // used to be getNodes(node,dir).size() but made it faster

template<size_t span>
void GraphUnitigsTemplate<span>::degree (const NodeGU& node, size_t &in, size_t &out) const  {  countNeighbors(node, in, out);  } 


template<size_t span>
int GraphUnitigsTemplate<span>::simplePathAvance (const NodeGU& node, Direction dir) const
{
    EdgeGU output;  return simplePathAvance (node, dir, output);
}

template<size_t span>
int GraphUnitigsTemplate<span>::simplePathAvance (const NodeGU& node, Direction dir, EdgeGU& output) const
{
    std::cout << "GraphU simplePathAvance called, not allowed in GraphUnitigs. try simplePathLongest_avance but it won't be edge-by-edge." << std::endl; exit(1);
    // shouldn't be called!
    // instead, do a more high level function
    return 0;
}


template<size_t span>
void GraphUnitigsTemplate<span>::simplify(unsigned int nbCores, bool verbose)
{
        Simplifications<GraphUnitigsTemplate<span>,NodeGU,EdgeGU> 
            graphSimplifications(*this, nbCores, verbose);
        graphSimplifications.simplify();
}


/*
 *
 * functions that aren't or shouldn't be implemented in GraphUnitigs
 *
 *
 */

template<size_t span> 
GraphIterator<NodeGU> GraphUnitigsTemplate<span>::getSimpleNodeIterator (const NodeGU& node, Direction dir) const
{
    std::cout << "getSimpleNodeIterator called in GraphU, not implemented" << std::endl; 
    exit(1);
    return GraphIterator<NodeGU> ();
}

template<size_t span> 
GraphIterator<EdgeGU> GraphUnitigsTemplate<span>::getSimpleEdgeIterator (const NodeGU& node, Direction dir) const
{
    std::cout << "getSimpleEdgeIterator called in GraphU, not implemented" << std::endl; 
    exit(1);
    return GraphIterator<EdgeGU>();
}

template<size_t span> 
int GraphUnitigsTemplate<span>::queryAbundance (const NodeGU& node) const
{
    std::cout << "queryAbundance called for a single node" << std::endl; 
    exit(1);
    return 0;
}

template<size_t span>
int GraphUnitigsTemplate<span>::queryNodeState (const NodeGU& node) const 
{
    std::cout << "queryNodeState called in GraphU" << std::endl; 
    std::cout << "not implemented" << std::endl; exit(1);
    return 0;
}

template<size_t span>
void GraphUnitigsTemplate<span>::setNodeState (const NodeGU& node, int state) const 
{
    //boost::apply_visitor (setNodeState_visitor<span>(node, state),  *(span*)BaseGraph::_variant);
    std::cout << "GraphUnitigs::setNodeState() not implemented" << std::endl; exit(1);
}

template<size_t span>
void GraphUnitigsTemplate<span>::resetNodeState() const
{
    //boost::apply_visitor (resetNodeState_visitor<span>(),  *(span*)BaseGraph::_variant);
    std::cout << "GraphUnitigs::resetNodeState() not implemented" << std::endl; exit(1);
}

template<size_t span>
void GraphUnitigsTemplate<span>::disableNodeState() const
{
    //boost::apply_visitor (disableNodeState_visitor<span>(),  *(span*)BaseGraph::_variant);
    std::cout << "GraphUnitigs::disableNodeState() not implemented" << std::endl;
}




/*
 *
 *
 * misc
 *
 *
 */

// instantiation
// uses Node and Edge as defined in Graph.hpp (legacy GATB compatibility, when Graph was not templated)
// so.. update;. this instantiation cannot be used by programs, but rather, it's to be used by TemplateSpecialization I think. So it could maybe be removed, and same for those in Graph.hpp
#ifdef WITH_LAMBDA_EXPRESSION //  requires C++11
template <size_t span>
using GraphUnitigs = GraphUnitigsTemplate<span>; 
#endif

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif