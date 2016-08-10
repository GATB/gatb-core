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

template<size_t span, typename Node_in>
unsigned long getNodeIndex (const GraphData<span>& data, Node_in& node);

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
    string unitigs_filename = "dummy.unitigs.fa"; // because there's already a bank, we don't know its name maybe?
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
    
    bool force_loading_unitigs = false;

    if (!checkState(STATE_BCALM2_DONE) && (!force_loading_unitigs /* for debug, if unitigs are made but the h5 didn't register it, stupid h5*/))
    {
        int nb_threads =
            props->getInt(STR_NB_CORES);

        size_t  kmerSize = BaseGraph::getKmerSize();
        if (kmerSize != (unsigned int)props->getInt(STR_KMER_SIZE))
            std::cout << "kmer discrepancy: should i take " << kmerSize << " or " << props->getInt(STR_KMER_SIZE) << std::endl;

        UnitigsConstructionAlgorithm<span> unitigs_algo(BaseGraph::getStorage(), unitigs_filename, nb_threads, props);

        BaseGraph::executeAlgorithm(unitigs_algo, &BaseGraph::getStorage(), props, BaseGraph::_info);
        
        setState(STATE_BCALM2_DONE);
    }

    /** We save the state at storage root level. */
    BaseGraph::getGroup().setProperty ("state",     Stringify::format("%d", BaseGraph::_state));
}


template<size_t span>
void GraphUnitigsTemplate<span>::load_unitigs(string unitigs_filename)
{
    BankFasta inputBank (unitigs_filename);
    //bank::IBank* inputBank = Bank::open (unitigs_filename);
    //LOCAL (inputBank);
    //ProgressIterator<bank::Sequence> itSeq (*inputBank, "loading unitigs");
    BankFasta::Iterator itSeq (inputBank);

    unsigned int kmerSize = BaseGraph::_kmerSize;

    uint32_t utig_counter = 0;
    uint64_t nb_utigs_nucl = 0;
    uint64_t nb_utigs_nucl_mem = 0;
    for (itSeq.first(); !itSeq.isDone(); itSeq.next()) // could be done in parallel, maybe, if we used many unordered_map's with a hash on the query kmer (TODO opt)
    {
        const string& seq = itSeq->toString();
        const string& comment = itSeq->getComment();
        float mean_abundance = atof(comment.substr(comment.find("MA=")+3).c_str());
    
        typename Model::Kmer kmerBegin = modelK->codeSeed(seq.substr(0, kmerSize).c_str(), Data::ASCII);
        typename Model::Kmer kmerEnd = modelK->codeSeed(seq.substr(seq.size() - kmerSize, kmerSize).c_str(), Data::ASCII);

        if (seq.size() > kmerSize)
        {
            ExtremityInfo eBegin(utig_counter, false, false, UNITIG_BEGIN);
            ExtremityInfo eEnd(utig_counter, false, false, UNITIG_END);
            utigs_map[kmerBegin.value()] = eBegin.pack();
            utigs_map[kmerEnd.value()] = eEnd.pack();
        }
        else // special case: single-kmer unitigs
        {
            ExtremityInfo eBegin(utig_counter, false, false, UNITIG_BOTH);
            utigs_map[kmerBegin.value()] = eBegin.pack();
        }
        unitigs.push_back(seq);
        unitigs_mean_abundance.push_back(mean_abundance);
        utig_counter++;
        nb_utigs_nucl += seq.size();
        nb_utigs_nucl_mem += seq.capacity();
    }

    unitigs_traversed.resize(0);
    unitigs_traversed.resize(utigs_map.size(), false); // resize "traversed" bitvector, setting it to zero as well


    //std::cout << "after load_unitigs utigs map size " << utigs_map.size() << std::endl;
//    utigs_map.rehash(0); // doesn't seem to help
    
    size_t count = 0;
    for (unsigned i = 0; i < utigs_map.bucket_count(); ++i) {
        size_t bucket_size = utigs_map.bucket_size(i);
        if (bucket_size == 0) {
            count++;
        }
        else {
            count += bucket_size;
        }
    }
    uint64_t utigs_map_capacity = count; // am not sure if that's really it. do values get allocated in an empty bucket?
    // anyway my tests show that i still underestimate the memory taken by the unordered_map
        

    // an estimation of memory usage
    uint64_t nb_kmers = utigs_map.size();
    uint64_t mem_vec = (unitigs.capacity() * sizeof(string) + nb_utigs_nucl_mem);
    std::cout << "Memory usage:" << std::endl;
    std::cout <<  "   " << (sizeof(Type) * utigs_map_capacity) / 1024 / 1024 << " MB keys in unitigs dict (Type = " << to_string(sizeof(Type)) << " bytes)" << std::endl;
    std::cout <<  "   " << (sizeof(uint32_t) * utigs_map_capacity) / 1024 / 1024 << " MB values in unitigs dict (values = " << to_string(sizeof(uint32_t)) << " bytes)" << std::endl;
    std::cout <<  "   " <<  mem_vec /1024 /1024 << " MB unitigs nucleotides" << std::endl;
    std::cout <<  "   " <<  (nb_kmers*4) / 1024 / 1024 << " MB unitigs abundances" << std::endl;
    std::cout <<  "   " <<  (utigs_map.size()/8) / 1024 / 1024 << " MB visited bitvector (hopefully 1 bit/elt)" << std::endl;
    std::cout <<  "Estimated total: " <<  (nb_kmers*(4*1.0/8.0) + utigs_map_capacity * ( sizeof(Type)  + sizeof(uint32_t)) + mem_vec) / 1024 / 1024 << " MB" << std::endl;

    std::cout << "utigs_map size calculated from buckets " << count << " vs size() " << utigs_map.size() << std::endl;
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
    
        utigs_map = graph.utigs_map;
        unitigs = graph.unitigs;
        unitigs_mean_abundance = graph.unitigs_mean_abundance;
        unitigs_traversed = graph.unitigs_traversed;
        modelK = graph.modelK;
        modelKdirect = graph.modelKdirect;
        
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
    
        utigs_map = std::move(graph.utigs_map);
        unitigs = std::move(graph.unitigs);
        unitigs_mean_abundance = std::move(graph.unitigs_mean_abundance);
        unitigs_traversed = std::move(graph.unitigs_traversed);
        modelK = std::move(graph.modelK);
        modelKdirect = std::move(graph.modelKdirect);
        
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

template<size_t span>
GraphVector<EdgeFast<span>> GraphUnitigsTemplate<span>::getEdges (NodeFast<span> source, Direction direction)  const
{
    bool debug = false;

    if (debug)
    std::cout << "graphU getEdges called" << std::endl;

    /* actually, nothing weird with having 1-in, 1-out nodes in GraphUnitigs: consider that example:
     *        --------
     *                 \
     *                  v
     *  -------[node] -> -----------
     *
     * [node] is 1-in and 1-out yet compactions were fine, it's just that there is in-branching in the following node.
     */
    //if (indegree(source) == 1 && outdegree(source) == 1)
 
    GraphVector<EdgeFast<span>> res;

    /*
    auto it = utigs_map.find(source.kmer);
    if (it == utigs_map.end())
    {   std::cout << std::endl << " source not found in utigs_map: " <<  BaseGraph::toString(source) << std::endl;
        throw system::Exception ("getEdges failure (see printed error above)");
    }*/

    const ExtremityInfo e(utigs_map.at(source.kmer)); // I used to have it->second here but for debugging, it's actually good that at() throws an exception

    bool same_orientation = node_in_same_orientation_as_in_unitig(source, e);
    res.resize(0);
            if (debug)
    std::cout << "source: " << BaseGraph::toString(source) << std::endl << "unitig: " << unitigs[e.unitig] << std::endl << "e: " << e.toString() << " same orientation " << same_orientation << std::endl;


    unsigned int kmerSize = BaseGraph::_kmerSize;
    const std::string& seq = unitigs[e.unitig];
    if (seq.size() > kmerSize)
    {
        // unitig: [kmer]-------
        if (same_orientation && (direction & DIR_OUTCOMING) && (e.pos & UNITIG_BEGIN)) 
        {
            Node dest = BaseGraph::buildNode(seq.substr(1, kmerSize).c_str());
            kmer::Nucleotide nt = BaseGraph::getNT(dest, kmerSize-1); 
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.kmer, source.strand, dest.kmer, dest.strand, nt, DIR_OUTCOMING);
            if (debug)
            std::cout << "found success of [kmer]---" << std::endl;
        }

        // unitig: [kmer rc]-------
        if ((!same_orientation) && (direction & DIR_INCOMING) && (e.pos & UNITIG_BEGIN)) 
        {
            Node dest = BaseGraph::buildNode(seq.substr(1, kmerSize).c_str());
            kmer::Nucleotide nt = BaseGraph::getNT(dest, kmerSize-1); 
            nt = reverse(nt); // TODO check it. i'm not sure about setting that nt. Minia doesn't use nt's so it might be wrong and i wouldn't know it
            dest = BaseGraph::reverse(dest);
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.kmer, source.strand, dest.kmer, dest.strand, nt, DIR_INCOMING);
            if (debug)
            std::cout << "found success of [kmer rc]---" << std::endl;
        }

        // unitig: ----------[kmer]
        if ((same_orientation) && (direction & DIR_INCOMING) && (e.pos & UNITIG_END))
        {
            Node dest = BaseGraph::buildNode(seq.substr(seq.size() - kmerSize - 1, kmerSize).c_str());
            kmer::Nucleotide nt = BaseGraph::getNT(dest, 0); 
            res.resize(res.size()+1);
            res[res.size()-1].set ( source.kmer, source.strand, dest.kmer, dest.strand, nt, DIR_INCOMING);
            if (debug)
            std::cout << "found predec of --------[kmer]" << std::endl;
        }

        // unitig: ----------[kmer rc]
        if ((!same_orientation) && (direction & DIR_OUTCOMING) && (e.pos & UNITIG_END))
        {
            Node dest = BaseGraph::buildNode(seq.substr(seq.size() - kmerSize - 1, kmerSize).c_str());
            kmer::Nucleotide nt = BaseGraph::getNT(dest, 0); 
            nt = reverse(nt); 
            res.resize(res.size()+1);
            dest = BaseGraph::reverse(dest);
            res[0].set ( source.kmer, source.strand, dest.kmer, dest.strand, nt, DIR_OUTCOMING);
            if (debug)
            std::cout << "found predec of --------[kmer rc]" << std::endl;
        }
    }
    
    // otherwise, that extremity kmer has neighbors at are also extremities.
    // so, mutate to get all 4 outneighrs, and test for their existence in the utigs_map

    if (debug)
        std::cout << "[out-of-unitig getEdges] for " << BaseGraph::toString(source) << " dir " << direction << " e: ["  << e.toString() << "]" << std::endl;
    
    bool incoming = false;

    auto functor = [&](const Type &neighbor){
        Type true_neighbor = neighbor;

        Type norm_neighbor = modelKdirect->reverse(true_neighbor);
        kmer::Strand strand = kmer::STRAND_REVCOMP;
        bool rc = true;
        if (true_neighbor < norm_neighbor) 
        {   
            rc = false;
            norm_neighbor = true_neighbor;
            strand = kmer::STRAND_FORWARD;
        }
        if (utigs_map.find(norm_neighbor) != utigs_map.end()) 
        {
            const ExtremityInfo e(utigs_map.at(norm_neighbor));
            if (e.deleted) 
            {
                if (debug)
                    std::cout << "found deleted neighbor unitig (kmer: " << modelKdirect->toString(norm_neighbor) << ")" << std::endl;
                return;
            }

            res.resize(res.size()+1);
            NodeFast<span> dest (norm_neighbor, strand);
            kmer::Nucleotide nt;
            if (!incoming)
                nt = BaseGraph::getNT(dest, 0); 
            else
                nt = BaseGraph::getNT(dest, kmerSize-1); 
            if (rc)
                nt = reverse(nt);
            
            if (debug)
            std::cout << "found kmer " << BaseGraph::toString(dest) << " (kmer: " << modelKdirect->toString(norm_neighbor) << " strand: " << dest.strand << ") dir " << toString(incoming?DIR_INCOMING:DIR_OUTCOMING) << " nt " << ascii(nt)  << std::endl;
            res[res.size()-1].set ( source.kmer, source.strand, dest.kmer, dest.strand, nt, (incoming?DIR_INCOMING:DIR_OUTCOMING));
        }
    }; 

    Type oriented_kmer = source.kmer;
    if (source.strand == kmer::STRAND_REVCOMP)
        oriented_kmer = modelKdirect->reverse(source.kmer);

    if (direction & DIR_OUTCOMING && ((same_orientation && (e.pos & UNITIG_END)) || ( !same_orientation && (e.pos & UNITIG_BEGIN)) ))
    {
        //modelKdirect->iterateOutgoingNeighbors(oriented_kmer, functor); // that's what i wanted to do, but it doesn't work, that function converts to canonical kmer anyhow
        /** We compute the 4 possible neighbors. */
        for (size_t nt=0; nt<4; nt++)
        {
            Type next1 = (((oriented_kmer) * 4 )  + nt) & modelKdirect->getKmerMax();
            functor(next1);
        }
    }
    if (direction & DIR_INCOMING && ((same_orientation && (e.pos & UNITIG_BEGIN)) || ( !same_orientation && (e.pos & UNITIG_END)) ))
    {
        incoming = true;

        Type rev = modelKdirect->reverse(oriented_kmer);
        for (size_t nt=0; nt<4; nt++)
        {
            Type next1 = (((rev) * 4 )  + nt) & modelKdirect->getKmerMax();
            functor(modelKdirect->reverse(next1));
        }
    }

   
    return res;
}

/* this function isn't the most efficient, but then again, Minia doesn't use it */
// I think I only implemented it so that it passes simple tests (TestDebruijnUnitigs.cpp)
template<size_t span>
GraphVector<NodeFast<span>> GraphUnitigsTemplate<span>::getNodes (NodeFast<span> &source, Direction direction)  const
{
    GraphVector<NodeFast<span>> nodes;
    GraphVector<EdgeFast<span>> edges = getEdges (source, direction);
    nodes.resize(edges.size());
    for (unsigned int i = 0; i < edges.size(); i++)
    {
        nodes[i] = edges[i].to;
    }
    return nodes;
}

// can be optimized surely, TODO opt
template<size_t span>
bool GraphUnitigsTemplate<span>::node_in_same_orientation_as_in_unitig(const Node& node, const ExtremityInfo& e) const
{
  int kmerSize = BaseGraph::_kmerSize;
  if (e.pos == UNITIG_BEGIN)
    return BaseGraph::toString(node) == unitigs[e.unitig].substr(0,kmerSize);
  else
    return BaseGraph::toString(node) == unitigs[e.unitig].substr(unitigs[e.unitig].size() - kmerSize, kmerSize);
}

template<size_t span>
unsigned char GraphUnitigsTemplate<span>::countNeighbors (NodeFast<span> &source, Direction direction)  const
{
    // for the sake of no duplication and removing bugs, i'm de-optimizing this function for now.
    GraphVector<EdgeFast<span>> edges = getEdges(source, direction);
    return edges.size();
}

template<size_t span>
void GraphUnitigsTemplate<span>::countNeighbors (NodeFast<span> &source, size_t &in, size_t &out)  const
{
    std::cout << "GraphU countNeighbors source,in,out not implememented" << std::endl;exit(1);
}

template<size_t span>
NodeFast<span> GraphUnitigsTemplate<span>::getNode (NodeFast<span>& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{
    std::cout << "GraphU getNode source,dir,nt,exists  not implememented" << std::endl;exit(1);
    return Node();
}

template<size_t span>
GraphIterator<NodeFast<span>> GraphUnitigsTemplate<span>::getNodes () const
{
    //std::cout << "graphUnitigs getNodes called" << std::endl;
    
    /* emulates iteration of nodes à la original GATB Graph */
    /* except that here, we only iterate the extremities of unitigs */
    class NodeIterator : public tools::dp::ISmartIterator<NodeFast<span>>
    {
        public:
            NodeIterator (const NodeMap& utigs_map) 
                :  utigs_map(utigs_map), _rank(0), _isDone(true)   {  
                    this->_item->strand = STRAND_FORWARD;  // iterated nodes are always in forward strand.
                    _nbItems = utigs_map.size();
                }

            ~NodeIterator ()  {  }

            u_int64_t rank () const { return _rank; }

            void update_item()
            {
                this->_rank ++;
                this->_item->kmer      = it->first;
                this->_item->abundance = 0; // returning per-node abundance in GraphUnitigs isn't supported
                const ExtremityInfo e(it->second);
                this->_item->mphfIndex = e.unitig; // rank of a node is GraphUnitigs is just its unitig index, _by convention_.
                this->_item->iterationRank = this->_rank;
            }

            /** \copydoc  Iterator::first */
            void first()
            {
                it = utigs_map.begin();
                _rank   = 0;
                _isDone = it == utigs_map.end();

                if (!_isDone)
                    update_item();
            }

            /** \copydoc  Iterator::next */
            void next()
            {
                it++;
                _isDone = it == utigs_map.end();
                if (!_isDone)
                    update_item();
            }

            /** \copydoc  Iterator::isDone */
            bool isDone() { return _isDone;  }

            /** \copydoc  Iterator::item */
            NodeFast<span>& item ()  {  return *(this->_item);  }

            /** is this ever called?*/
            void setItem (NodeFast<span>& i)
            {
                /** We set the node item to be set for the current iterator. */
                this->_item = &i;
                this->_item->strand = STRAND_FORWARD;
            }

            /** */
            u_int64_t size () const { return _nbItems; }

        private:
            const NodeMap& utigs_map;
            typename NodeMap::const_iterator it;

            u_int64_t _rank;
            bool      _isDone;
            u_int64_t _nbItems;
    };

    return new NodeIterator(utigs_map);
}

template<size_t span> 
GraphIterator<NodeFast<span>> GraphUnitigsTemplate<span>::getSimpleNodeIterator (NodeFast<span>& node, Direction dir) const
{
    std::cout << "getSimpleNodeIterator called in GraphU, not implemented" << std::endl; 
    exit(1);
    return GraphIterator<NodeFast<span>> ();
}

template<size_t span> 
GraphIterator<EdgeFast<span>> GraphUnitigsTemplate<span>::getSimpleEdgeIterator (NodeFast<span>& node, Direction dir) const
{
    std::cout << "getSimpleEdgeIterator called in GraphU, not implemented" << std::endl; 
    exit(1);
    return GraphIterator<EdgeFast<span>>();
}


/** */
template<size_t span> 
int GraphUnitigsTemplate<span>::queryAbundance (NodeFast<span>& node) const
{
    std::cout << "queryAbundance called for a single node" << std::endl; 
    exit(1);
    return 0;
}

/** */
template<size_t span>
int GraphUnitigsTemplate<span>::queryNodeState (NodeFast<span>& node) const 
{
    std::cout << "queryNodeState called in GraphU" << std::endl; 
    std::cout << "not implemented" << std::endl; exit(1);
    return 0;
}


/** */
template<size_t span>
void GraphUnitigsTemplate<span>::setNodeState (NodeFast<span>& node, int state) const 
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

/* this can be useful for benchmarking purpose.
 * because contains() always does a MPHF query to check if the node is deleted or not
 * unless _nodestate is NULL.
 * thus, this functions sets _nodestate to NULL.
 * NOTE: irreversible!
 */
template<size_t span>
void GraphUnitigsTemplate<span>::disableNodeState() const
{
    //boost::apply_visitor (disableNodeState_visitor<span>(),  *(span*)BaseGraph::_variant);
    std::cout << "GraphUnitigs::disableNodeState() not implemented" << std::endl;
}

template<size_t span> 
bool GraphUnitigsTemplate<span>::isNodeDeleted(NodeFast<span>& node) const
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    return e.deleted;
}


// direct access to MPHF index

template<size_t span> 
unsigned long GraphUnitigsTemplate<span>::nodeMPHFIndex(NodeFast<span>& node) const 
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    return (((e.pos & UNITIG_END)) ? 1 : 0)* utigs_map.size() + e.unitig; 
}

/* warning: will delete the whole simple path */
template<size_t span>
void GraphUnitigsTemplate<span>::deleteNode (NodeFast<span>& node) 
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
bool GraphUnitigsTemplate<span>::contains (const NodeFast<span>& item) const
{
    std::cout << "contains() not implemeneted in GraphUnitigs" << std::endl; exit(1);
    return false;
}

template<size_t span>
bool GraphUnitigsTemplate<span>::contains (const typename Kmer<span>::Type& item) const
{
    std::cout << "contains() not implemeneted in GraphUnitigs" << std::endl; exit(1);
}



// high-level functions that used to be in Simplifications.cpp

template<size_t span>
bool GraphUnitigsTemplate<span>::
isLastNode(const NodeFast<span>& node, Direction dir) const
{
    if (utigs_map.find(node.kmer) == utigs_map.end()) { std::cout << "error: node not indexed " << BaseGraph::toString(node) << std::endl; exit(1); }

    const ExtremityInfo e(utigs_map.at(node.kmer));
    bool same_orientation = node_in_same_orientation_as_in_unitig(node, e);

    // cases where, following that unitig direction, we're already at the last node
    if ((same_orientation    && (e.pos & UNITIG_END) && dir == DIR_OUTCOMING) ||
        (same_orientation    && (e.pos & UNITIG_BEGIN) && dir == DIR_INCOMING) ||
        ((!same_orientation) && (e.pos & UNITIG_END) && dir == DIR_INCOMING) ||
        ((!same_orientation) && (e.pos & UNITIG_BEGIN) && dir == DIR_OUTCOMING))
        return true;
    return false;
}

template<size_t span>
bool GraphUnitigsTemplate<span>::
isFirstNode(const NodeFast<span>& node, Direction dir) const
{
    return !isLastNode(node,dir);
}

template<size_t span>
double GraphUnitigsTemplate<span>::
simplePathMeanAbundance     (const NodeFast<span>& node, Direction dir) 
{
    if (isLastNode(node,dir))
        return 0;

    float coverage = 0;
    int endDegree;
    int seqLength = 0;
    simplePathLongest_avance(node, dir, seqLength, endDegree, false /*markDuringTraversal*/, coverage);
    return coverage / (float)seqLength;
}

template<size_t span>
double GraphUnitigsTemplate<span>::
unitigMeanAbundance     (const NodeFast<span>& node) const
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    return unitigs_mean_abundance[e.unitig];
}

template<size_t span>
unsigned int GraphUnitigsTemplate<span>::
simplePathLength            (const NodeFast<span>& node, Direction dir) 
{
    float coverage = 0;
    int endDegree;
    int seqLength = 0;
    simplePathLongest_avance(node, dir, seqLength, endDegree, false /*markDuringTraversal*/, coverage);
    return seqLength;
}

template<size_t span>
unsigned int GraphUnitigsTemplate<span>::
unitigLength            (const NodeFast<span>& node, Direction dir) const
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    const std::string& seq = unitigs[e.unitig];
    bool same_orientation = node_in_same_orientation_as_in_unitig(node, e);

    int length;
    if (isLastNode(node,dir))
        length = 0;
    else
        length = unitigs[e.unitig].size() - BaseGraph::_kmerSize;
    return length;
}

template<size_t span>
NodeFast<span> GraphUnitigsTemplate<span>::
unitigLastNode          (const NodeFast<span>& node, Direction dir) const
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    const std::string& seq = unitigs[e.unitig];

    //std::cout << "lastnode, same orientation " << same_orientation << " e " <<  e.toString() << " dir " << dir  << std::endl;
    
    if (isLastNode(node,dir))
        return node;

    int kmerSize = BaseGraph::_kmerSize;
    NodeFast<span> res;
    if ((e.pos & UNITIG_BEGIN))
        res = BaseGraph::buildNode(seq.substr(seq.size() - kmerSize, kmerSize).c_str());
    else
        res = BaseGraph::buildNode(seq.substr(0, kmerSize).c_str());
    bool same_orientation = node_in_same_orientation_as_in_unitig(node, e);
    if (!same_orientation)
        res = BaseGraph::reverse(res);
    return res;
}

template<size_t span>
NodeFast<span> GraphUnitigsTemplate<span>::
simplePathLastNode          (const NodeFast<span>& node, Direction dir) 
{
    std::vector<NodeFast<span>> nodesList;
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
unitigDelete (NodeFast<span>& node) 
{
    ExtremityInfo e(utigs_map[node.kmer]);
    e.deleted = true;
    utigs_map[node.kmer] = e.pack();
    
    
    // make sure to delete other part also
    const std::string& seq = unitigs[e.unitig];
    int kmerSize = BaseGraph::_kmerSize;
    NodeFast<span> second;
    if ((e.pos & UNITIG_BEGIN)) // TODO might make sense to replace all this by get_other_end(Node,dir) with a check that we're not in the wrong direction. but not needed in minia so i didn't do it.
        second = BaseGraph::buildNode(seq.substr(seq.size() - kmerSize, kmerSize).c_str());
    else
        second = BaseGraph::buildNode(seq.substr(0, kmerSize).c_str());

    ExtremityInfo e2(utigs_map[second.kmer]);
    e2.deleted = true;
    utigs_map[second.kmer] = e2.pack();

    //std::cout << "GraphU deleted simple path " << seq << std::endl; 
}


// wrapper that only puts the node into nodesDeleter
template<size_t span>
void GraphUnitigsTemplate<span>::
unitigDelete (NodeFast<span>& node, Direction dir, NodesDeleter<NodeFast<span>, EdgeFast<span>, GraphUnitigsTemplate<span>>& nodesDeleter) 
{
    nodesDeleter.onlyListMethod = true; // a bit inefficient to always tell the deleter to be in that mode, but so be it for now. just 1 instruction, won't hurt.
    nodesDeleter.markToDelete(node);
}

template<size_t span>
void GraphUnitigsTemplate<span>::
simplePathDelete (NodeFast<span>& node, Direction dir, NodesDeleter<NodeFast<span>, EdgeFast<span>, GraphUnitigsTemplate<span>>& nodesDeleter) 
{
    std::vector<NodeFast<span>> nodesList;
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
unitigSequence (const NodeFast<span>& node, bool& isolatedLeft, bool& isolatedRight) const
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    const string& seq = unitigs[e.unitig];

    //std::cout << " seq " << seq << " node " << BaseGraph::toString(node) << std::endl;
    int kmerSize = BaseGraph::_kmerSize;
    Node right = BaseGraph::buildNode(seq.substr(seq.size() - kmerSize, kmerSize).c_str());
    Node left = BaseGraph::buildNode(seq.substr(0, kmerSize).c_str());

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
simplePathLongest_avance(const NodeFast<span>& node, Direction dir, int& seqLength, int& endDegree, bool markDuringTraversal, float& coverage, string* seq, std::vector<NodeFast<span>> *nodesList) 
{
    bool debug = false;
    if (debug)
        std::cout << "simplePathLongest_avance called, node " << BaseGraph::toString(node) << " dir " << dir << std::endl;

    seqLength = 0;
    int kmerSize = BaseGraph::_kmerSize;
    NodeFast<span> cur_node = node;
    if (isFirstNode(cur_node, dir))
    {

        // first node in unitig may have in-branching, it's fine for a simple path traversal. we'll just go to last node and record the sequence of that unitig
        const ExtremityInfo e(utigs_map[node.kmer]);
        string new_seq = unitigs[e.unitig];

        bool same_orientation = node_in_same_orientation_as_in_unitig(node, e);
        if (!same_orientation)
            new_seq = revcomp(new_seq);

        if (seq != nullptr)
        {
            if (dir == DIR_OUTCOMING)
                *seq += new_seq.substr(kmerSize-1);
            else
                *seq = new_seq.substr(0,new_seq.size()-(kmerSize-1)) + *seq;
        }

        // length is all the kmers of that unitig, except first one
        seqLength += new_seq.size() - kmerSize;
        // for coverage: it includes the abundance of the first kmer, because there's no way to exclude it, hence the +1
        coverage += unitigMeanAbundance(cur_node) * (new_seq.size() - kmerSize + 1);
       
        if (debug)
            std::cout << "simplePathLongest_avance was at a first node = " << BaseGraph::toString(cur_node) << " strand " << cur_node.strand << " so traversed unitig of length " << new_seq.size() << std::endl;

        cur_node = unitigLastNode(node,dir);       
        
        if (debug)
            std::cout << "simplePathLongest_avance now at last node = " << BaseGraph::toString(cur_node) << " strand " << cur_node.strand << std::endl;

        if (nodesList != nullptr)
            nodesList->push_back(cur_node);
    }

    if (markDuringTraversal) // no need to mark, that unitig should already be marked by minia, but doing it anyway just to be safe
        unitigMark(cur_node);

    while (true)
    { // invariant here: cur_node is the last node of a unitig

        assert(isLastNode(cur_node,dir));

        GraphVector<EdgeFast<span>> neighbors = this->neighborsEdge (cur_node, dir);
        endDegree = neighbors.size();
        
        /** We check we have no outbranching. */
        if (endDegree != 1)
        {
            if (debug)
                std:: cout << "simplePathLongest_avance stopped because # neighbor of node " << BaseGraph::toString(cur_node) << ": " << neighbors.size() << std::endl;
            return;
        }
      
        const ExtremityInfo e(utigs_map[neighbors[0].to.kmer]);
        string new_seq = unitigs[e.unitig];

        bool same_orientation = node_in_same_orientation_as_in_unitig(neighbors[0].to, e);

        if (debug)
            std::cout << "simplePathLongest_avance continues now at a last node = " << BaseGraph::toString(cur_node) << " strand " << cur_node.strand << " neighbor.to " << BaseGraph::toString(neighbors[0].to) << " strand " << neighbors[0].to.strand  << " new seq length: " << new_seq.size() << std::endl;


        if (!same_orientation)
            new_seq = revcomp(new_seq);

        // some sanity checks      
        if (e.pos != UNITIG_BOTH) 
        {
            // FIXME: there is a bcalm bug. see strange_seq.fa. i'll send it to antoine, but meanwhile, i'm coding a workaround
            if (0)
            {
            if (dir==DIR_INCOMING )
                assert( (e.pos == UNITIG_END   && same_orientation) || (e.pos == UNITIG_BEGIN && (!same_orientation)));
            else
                assert( (e.pos == UNITIG_BEGIN && same_orientation) || (e.pos == UNITIG_END   && (!same_orientation)));
            }
            if (dir==DIR_INCOMING )
            {
                if (!( (e.pos == UNITIG_END   && same_orientation) || (e.pos == UNITIG_BEGIN && (!same_orientation))))
                {
                    unitigDelete(neighbors[0].to);
                    return;
                }
            }
            else
            {
                if (!( (e.pos == UNITIG_BEGIN && same_orientation) || (e.pos == UNITIG_END   && (!same_orientation))))
                {
                    unitigDelete(neighbors[0].to);
                    return;
                }
            }

        }
       
        GraphVector<EdgeFast<span>> in_neighbors_vec = this->neighborsEdge (neighbors[0].to, reverse(dir));
        int in_neighbors = in_neighbors_vec.size();

        assert(in_neighbors >= 1);
        /** We check we have no in-branching. */
        if (in_neighbors > 1) 
        {
            if (debug)
                std:: cout << "simplePathLongest_avance stopped at " << BaseGraph::toString(cur_node) << " because of in-branching " << in_neighbors << std::endl;
            return;
        } 

        NodeFast<span> last_node = unitigLastNode(neighbors[0].to, dir);
        cur_node = last_node;

        if (nodesList != nullptr)
            nodesList->push_back(cur_node);

        // append the sequence (except the overlap part, of length k-1.
        if (seq != nullptr)
        {
            if (dir == DIR_OUTCOMING)
                *seq += new_seq.substr(kmerSize-1);
            else
                *seq = new_seq.substr(0,new_seq.size()-(kmerSize-1)) + *seq;
        }

        seqLength += new_seq.size() - (kmerSize-1);
        coverage += unitigMeanAbundance(cur_node) * (new_seq.size() - kmerSize + 1); // here too, coverage is computed according to whole unitig

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
simplePathBothDirections(const NodeFast<span>& node, bool& isolatedLeft, bool& isolatedRight, bool markDuringTraversal, float &coverage) 
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    string seq = unitigs[e.unitig];
    
    int kmerSize = BaseGraph::_kmerSize;
    float midTotalCoverage = unitigMeanAbundance(node) * (seq.size() - kmerSize + 1);

    Node left = BaseGraph::buildNode(seq.substr(0, kmerSize).c_str());
    Node right = BaseGraph::buildNode(seq.substr(seq.size() - kmerSize, kmerSize).c_str());

    //std::cout << "starting seq " << seq << " (from node " << BaseGraph::toString(node) << ") left: " << BaseGraph::toString(left) << " right: " << BaseGraph::toString(right) << std::endl;

    if (markDuringTraversal)
        unitigMark(left);

    string seqRight, seqLeft;
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
unitigMark            (const Node& node) 
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    unitigs_traversed[e.unitig] = true;
} 

template<size_t span>
bool GraphUnitigsTemplate<span>::
unitigIsMarked        (Node& node) const 
{
    const ExtremityInfo e(utigs_map.at(node.kmer));
    return unitigs_traversed[e.unitig];
}

/* when debugging goes wrong;. print the whole graph. */
template<size_t span>
void GraphUnitigsTemplate<span>::
debugPrintAllUnitigs() const
{
    std::cout << "Debug: printing all graph unitigs and status" << std::endl;
    for (auto utig : utigs_map)
    {
        const ExtremityInfo e(utig.second);
        std::cout << "unitig " << e.unitig << " extremity " << BaseGraph::toString(utig.first) << " " << e.toString() << " unitig " << unitigs[e.unitig] << " (length: " << unitigs[e.unitig].size() << ")" << " links: ";

        // printing neighbors
        string nodeStr = modelKdirect->toString(utig.first);
        Node cur_node = BaseGraph::buildNode(nodeStr.c_str());

        for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
        { 
            if (dir == DIR_OUTCOMING)
                std::cout << "out: ";
            else
                std::cout << "in: ";
            GraphVector<EdgeFast<span>> neighbors = this->neighborsEdge (cur_node, dir);
            for (size_t i = 0; i < neighbors.size(); i++)
            {
                if ( utigs_map.find(neighbors[i].to.kmer) == utigs_map.end())
                {
                    std::cout << "[in-unitig]";
                    break; // in that direction there's nothing
                }

                const uint32_t packed = utigs_map.at(neighbors[i].to.kmer);
                const ExtremityInfo en(packed);
#if 0
                string nodeStr;
                if (en.pos == UNITIG_BEGIN)
                nodeStr =  unitigs[en.unitig].substr(0,kmerSize);
                else
                nodeStr = unitigs[en.unitig].substr(unitigs[en.unitig].size() - kmerSize, kmerSize);
#endif
                string pos = (en.pos == UNITIG_BEGIN)?"beg":((en.pos == UNITIG_BOTH)?"kmer":"end");
                std::cout << en.unitig << "(" << pos << ") ";
            }
            std::cout <<  "    ";
        }
            std::cout << std::endl;
    }
    std::cout << "Done printing graph unitigs" << std::endl;
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
 * no changes needed in GraphUnitigs
 *
 */


template<size_t span>
bool GraphUnitigsTemplate<span>::isBranching (NodeFast<span>& node) const
/** REMARKS : it's a bit of a misnomer, but in GATB language, Branching means "anything not simple". I.e. simple means outdegree=indegree=1 and branching is every other degree combination, including 0's */
{
    size_t in, out;
    degree(node, in, out); 
    return (! (in==1 && out==1));
}

template<size_t span>
bool GraphUnitigsTemplate<span>::isSimple (Edge& edge) const
{    return this->outdegree(edge.from)==1  &&  this->indegree(edge.to)==1;}

template<size_t span>
size_t GraphUnitigsTemplate<span>::indegree (NodeFast<span>& node) const  {  return degree(node, DIR_INCOMING);   }

template<size_t span>
size_t GraphUnitigsTemplate<span>::outdegree (NodeFast<span>& node) const  {  return degree(node, DIR_OUTCOMING);  }

template<size_t span>
size_t GraphUnitigsTemplate<span>::degree (NodeFast<span>& node, Direction dir) const  {  return countNeighbors(node, dir);  } // used to be getNodes(node,dir).size() but made it faster

template<size_t span>
void GraphUnitigsTemplate<span>::degree (NodeFast<span>& node, size_t &in, size_t &out) const  {  countNeighbors(node, in, out);  } 

template<size_t span>
int GraphUnitigsTemplate<span>::simplePathAvance (NodeFast<span>& node, Direction dir, kmer::Nucleotide& nt) const
{
    Edge edge;
    int res = simplePathAvance (node, dir, edge);
    nt = edge.nt;
    return res;
}

template<size_t span>
int GraphUnitigsTemplate<span>::simplePathAvance (NodeFast<span>& node, Direction dir) const
{
    Edge output;  return simplePathAvance (node, dir, output);
}

template<size_t span>
int GraphUnitigsTemplate<span>::simplePathAvance (NodeFast<span>& node, Direction dir, Edge& output) const
{
    std::cout << "GraphU simplePathAvance called, not allowed in GraphUnitigs. try simplePathLongest_avance but it won't be edge-by-edge." << std::endl; exit(1);
    // shouldn't be called!
    // instead, do a more high level function

    GraphVector<EdgeFast<span>> neighbors = this->neighborsEdge (node, dir);

    /** We check we have no outbranching. */
    if (neighbors.size() == 1)
    {
        /** We set the output result. */
        output = neighbors[0];

        /** We check whether the neighbor has an inbranching or not. */
        if (this->degree (neighbors[0].to, impl::reverse(dir)) > 1)  {  return -2;  }

        /** We have a simple node here. */
        return 1;
    }

    /** if this kmer has out-branching, don't extend it. */
    if (neighbors.size() > 1)  {  return -1;  }

    return 0;
}

// could be factorized with Graph.hpp somehow
template<size_t span>
GraphVector<EdgeFast<span>> GraphUnitigsTemplate<span>::getEdgeValues (const typename NodeFast<span>::Value& kmer) const
{
    NodeFast<span> source (kmer);

    GraphVector<EdgeFast<span>> v1 = getEdges(source,          DIR_OUTCOMING);
    GraphVector<EdgeFast<span>> v2 = getEdges(BaseGraph::reverse(source), DIR_OUTCOMING);
#if 0
    v1.insert (v1.end(), v2.begin(), v2.end());
#else
    size_t n1=v1.size(), n2=v2.size();
    v1.resize (n1+n2);
    for (size_t i=0; i<n2; i++)  { v1[i+n1] = v2[i];  }
#endif

    return v1;
}


template<size_t span>
void GraphUnitigsTemplate<span>::simplify(unsigned int nbCores, bool verbose)
{
        Simplifications<GraphUnitigsTemplate<span>,NodeFast<span>,EdgeFast<span>> 
            graphSimplifications(*this, nbCores, verbose);
        graphSimplifications.simplify();
}


/*
 *
 * functions that shouldn't be implemented in GraphUnitigs
 *
 *
 */



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
