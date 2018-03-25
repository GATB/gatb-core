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

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/BranchingAlgorithm.hpp>
#include <gatb/debruijn/api/IContainerNode.hpp>

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
#include <gatb/kmer/impl/BloomAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithmFactory.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/CountProcessor.hpp>
#include <gatb/kmer/impl/MPHFAlgorithm.hpp>
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

#ifndef _GATB_CORE_DEBRUIJN_IMPL_GRAPHCPP_
#define _GATB_CORE_DEBRUIJN_IMPL_GRAPHCPP_

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

template<size_t span, typename Node_in>
unsigned long getNodeIndex (const GraphData<span>& data, Node_in& node);

/** These two functions set a specific value to the given GraphDataVariant object,
 * according to the provided kmer size.
 *
 * This method should be called in each Graph constructor in order to be sure to have
 * its data variant attribute set with the correct type.
 *
 * Many methods of the Graph class use the boost::static_visitor mechanism in order to
 * retrieve the "correct" type pointed by the boost variant. So it is mandatory that
 * each Graph instance has its data variant correctly set thanks to this "setVariant"
 * function.
 */
template<size_t span> struct FunctorSetVariant
{
    // actually this doesn't seem to help
    //template<typename GraphDataVariant>
    void operator ()  (GraphDataVariant& data)  {  data = GraphData<span> ();  }
};

template<typename Node, typename Edge, typename GraphDataVariant>
struct setVariant_visitor : public boost::static_visitor<>    {
    setVariant_visitor () {}

    template<size_t span>  void operator() (GraphData<span>& data) const
    {
        data = GraphData<span>();
    }
};
 
template<typename Node, typename Edge, typename GraphDataVariant_t>
void GraphTemplate<Node, Edge, GraphDataVariant_t>::setVariant (void* data, size_t kmerSize, size_t integerPrecision)
{
	/** We may force the kmer size and not use the optimized KSIZE value. */
	if (integerPrecision > 0)  { kmerSize = integerPrecision*32 - 1; }

    if (typeid(GraphDataVariant_t) == typeid(GraphDataVariant))
        Integer::apply<FunctorSetVariant,GraphDataVariant&> (kmerSize, (GraphDataVariant&)(*((GraphDataVariant*)data)));
    else
        // boost::apply_visitor (setVariant_visitor<Node, Edge, GraphDataVariant_t>(),  (GraphDataVariant_t&)data); // it is strange that this line doesn't work
        boost::apply_visitor (setVariant_visitor<Node, Edge, GraphDataVariant_t>(),  *(GraphDataVariant_t*)data); // but this one does
        //*((GraphDataVariant_t*)data) = (GraphDataVariant_t) (GraphData<32>()); // this one as well, but only for a single span
}



/********************************************************************************/

template<size_t span>
struct Count2TypeAdaptor  {  typename Kmer<span>::Type& operator() (typename Kmer<span>::Count& c)  { return c.value; }  };

/* This visitor is used to configure a GraphDataVariant object (ie configure its attributes).
 * The information source used to configure the variant is a kmer size and a storage.
 *
 * If the storage is null, only the kmer model is set.
 *
 * If the storage is not null, it is likely coming from a previous graph building (dsk, debloom, ...); we can
 * therefore configure the variant with the items coming from this storage.
 */
template<typename Node, typename Edge, typename GraphDataVariant>
template <size_t span>  
void configure_visitor<Node,Edge,GraphDataVariant>::operator() (GraphData<span>& data) const 
{
    size_t   kmerSize = graph.getKmerSize();
    
    /** We create the kmer model. */
    data.setModel (new typename Kmer<span>::ModelCanonical (kmerSize));

    if (graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_CONFIGURATION_DONE)
    {
        /** We get the configuration group in the storage. */
        Group& configGroup = storage.getGroup("configuration");

        /** We read the XML file and update the global info. */
        stringstream ss; ss << configGroup.getProperty ("xml");

        Properties props; props.readXML (ss);
        graph.getInfo().add (1, props);
    }

    if (graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE)
    {
        /** We get the dsk group in the storage. */
        Group& dskGroup = storage.getGroup("dsk");

        /** We set the iterable for the solid kmers. */
        data.setSolid (& dskGroup.getPartition<typename Kmer<span>::Count> ("solid"));

        /** We read the XML file and update the global info. */
        stringstream ss; ss << dskGroup.getProperty ("xml");

        Properties props; props.readXML (ss);
        graph.getInfo().add (1, props);
    }

    if (graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BLOOM_DONE)
    {
        /** We set the container. */
        BloomAlgorithm<span> algo (storage);
        graph.getInfo().add (1, algo.getInfo());
    }

    if (graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_DEBLOOM_DONE)
    {
        /** We set the container. */
        DebloomAlgorithm<span> algo (storage);
        graph.getInfo().add (1, algo.getInfo());
        data.setContainer (algo.getContainerNode());
    }

    if (graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BRANCHING_DONE)
    {
        /** We set the branching container. */
        BranchingAlgorithm<span, Node, Edge, GraphTemplate<Node, Edge, GraphDataVariant>> algo (storage);
        graph.getInfo().add (1, algo.getInfo());
        data.setBranching (algo.getBranchingCollection());
    }

    if ((graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE) &&  (graph.getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE))
    {
        typedef typename Kmer<span>::Count Count;
        typedef typename Kmer<span>::Type  Type;

        /** We get the dsk group in the storage. */
        Group& dskGroup = storage.getGroup("dsk");

        /** We get the iterable for the solid counts and solid kmers. */
        Partition<Count>* solidCounts = & dskGroup.getPartition<Count> ("solid");
        Iterable<Type>*   solidKmers  = new IterableAdaptor<Count,Type,Count2TypeAdaptor<span> > (*solidCounts);

        MPHFAlgorithm<span> mphf_algo (
                dskGroup,
                "mphf",
                solidCounts,
                solidKmers,
                1,  // loading using 1 thread
                false  /* build=true, load=false */
                );

        data.setAbundance (mphf_algo.getAbundanceMap());
        data.setNodeState (mphf_algo.getNodeStateMap());
        data.setAdjacency (mphf_algo.getAdjacencyMap());
    }
}


/* used in build_visitor and build_visitor_postsolid */
/** Algorithm configuration. */
template<typename Node, typename Edge, typename GraphDataVariant_t>
void GraphTemplate<Node, Edge, GraphDataVariant_t>::
executeAlgorithm (Algorithm& algorithm, Storage* storage, IProperties* props, IProperties& info) 
{
    algorithm.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));

    algorithm.run ();

    info.add (1, algorithm.getInfo());
    info.add (1, algorithm.getSystemInfo());

    if (storage != 0)
    {
        /** We memorize information of the algorithm execution as a property of the corresponding group. */
        storage->getGroup(algorithm.getName()).setProperty("xml", string("\n") + algorithm.getInfo()->getXML());
    }
}


/********************************************************************************/

/* These two visitors are used to build a graph. In particular, the data variant of the graph will
 * be configured through this boost visitor.
 *
 * The skeleton of the graph building is the following:
 *
 *  build_visitor_solid:
 *  - conversion of the input reads into a binary format
 *  - kmers counting
 *
 *  build_visitor_postsolid
 *  - mphf
 *  - deblooming
 *  - branching nodes computation
 *
 *  The common source for storing results for these different parts is a Storage instance, configured
 *  at the beginning of this visitor. This common storage is then provided to the different parts.
 *
 *  The data variant itself is configured throughout the steps. After that, the graph has enough
 *  knowledge to provide all the services of the Graph API. In particular, the branching nodes computation
 *  is based on the Graph API methods, and is therefore executed after the data variant configuration.

 */

template<typename Node, typename Edge, typename GraphDataVariant>
template <size_t span>  
void build_visitor_solid<Node,Edge,GraphDataVariant>::operator() (GraphData<span>& data) const 
{
    /** Shortcuts. */
    typedef typename Kmer<span>::Count Count;

    LOCAL (bank);

    // TODO Erwan: do we even use those variables? if not, why put them here? I feel that it's weak to parse cmdline arguments in build_visitor
    // because graph._kmerSize is already defined at this point; and config._minim_size has minimizer size info. The rest we don't use.
    size_t kmerSize      = props->get(STR_KMER_SIZE)          ? props->getInt(STR_KMER_SIZE)           : 31;
    size_t compressLevel   = props->get(STR_COMPRESS_LEVEL)   ? props->getInt(STR_COMPRESS_LEVEL)      : 0;
    bool   configOnly      = props->get(STR_CONFIG_ONLY) != 0;

    string output = props->get(STR_URI_OUTPUT) ?
        props->getStr(STR_URI_OUTPUT)   :
        (props->getStr(STR_URI_OUTPUT_DIR) + "/" + system::impl::System::file().getBaseName (bank->getId()));

    /* create output dir if it doesn't exist */
    if(!System::file().doesExist(props->getStr(STR_URI_OUTPUT_DIR))){
        int ok = System::file().mkdir(props->getStr(STR_URI_OUTPUT_DIR), 0755);
        if(ok != 0){
            throw system::Exception ("Error: can't create output directory");
        }
    }

    DEBUG ((cout << "builGraph for bank '" << bank->getId() << "'"
                << " kmerSize=" << kmerSize
                << " output='" << output << "'"
                << endl
           ));

    /** We create the kmer model. */
    data.setModel (new typename Kmer<span>::ModelCanonical (kmerSize));

    /** We add library and host information. */
    graph.getInfo().add (1, & LibraryInfo::getInfo());
    graph.getInfo().add (1, & HostInfo::getInfo());

    /************************************************************/
    /*                       Storage creation                   */
    /************************************************************/

    Storage* mainStorage = StorageFactory(graph._storageMode).create (output, true, false); /* second arg true = delete if exists; we're recreating this hdf5 file */

    /** We create the storage object for the graph. */
    graph.setStorage (mainStorage);

    /** We may need a specific storage instance if we want to save the solid kmers in a separate file. */
    Storage* solidStorage = 0;
    if (props->get(STR_URI_SOLID_KMERS) != 0)
    {
        string solidsName = props->getStr(STR_URI_SOLID_KMERS);
        /** Historically (by convention), the file was deleted if its name is "none".
         * Now, since debloom may use the minimizer repartition function (stored in the solid file),
         * we must not delete the solid file. */
        bool autoDelete = false; // (solidsName == "none") || (solidsName == "null");
        solidStorage = StorageFactory(graph._storageMode).create (solidsName, true, autoDelete);
    }
    else
    {
        solidStorage = graph._storage;
    }
    LOCAL (solidStorage);

    /** We change the compression level for the storage.  Note that we need to do this before accessing the groups. */
    mainStorage->root(). setCompressLevel (compressLevel);
    solidStorage->root().setCompressLevel (compressLevel);

    /** We get the minimizers hash group in the storage object. */
    Group& minimizersGroup = (*mainStorage)("minimizers");

    /** We get the 'dsk' group in the storage object. */
    Group& dskGroup = (*solidStorage)("dsk");

    /************************************************************/
    /*                       Configuration                      */
    /************************************************************/
    DEBUG ((cout << "build_visitor : ConfigurationAlgorithm BEGIN\n"));

    ConfigurationAlgorithm<span> configAlgo (bank, props);
    configAlgo.getInput()->add (0, STR_STORAGE_TYPE, std::to_string(graph._storageMode) );
    graph.executeAlgorithm (configAlgo, & graph.getStorage(), props, graph._info);
    Configuration config = configAlgo.getConfiguration();
    graph.setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_CONFIGURATION_DONE);

    /* remember configuration details (e.g. number of passes, partitions). useful for bcalm. */                                                                     
    graph.getStorage().getGroup(configAlgo.getName()).setProperty("xml", string("\n") + configAlgo.getInfo()->getXML());        

    DEBUG ((cout << "build_visitor : ConfigurationAlgorithm END\n"));

    /** We may have to stop just after configuration. */
    if (configOnly)  { return; }

    /************************************************************/
    /*                  Minimizers repartition                  */
    /************************************************************/
    DEBUG ((cout << "build_visitor : RepartitorAlgorithm BEGIN\n"));

    RepartitorAlgorithm<span> repart (bank, 
            minimizersGroup, 
            config,
            props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
            );
    graph.executeAlgorithm (repart, 0, props, graph._info);

    DEBUG ((cout << "build_visitor : RepartitorAlgorithm END\n"));

    /************************************************************/
    /*                         Sorting count                    */
    /************************************************************/
    DEBUG ((cout << "build_visitor : SortingCountAlgorithm BEGIN\n"));

    /** We create a DSK instance and execute it. */
    SortingCountAlgorithm<span> sortingCount (
            bank,
            config,
            new Repartitor(minimizersGroup),
            SortingCountAlgorithm<span>::getDefaultProcessorVector (config, props, solidStorage, mainStorage),
            props
            );

    graph.executeAlgorithm (sortingCount, solidStorage, props, graph._info);
    graph.setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE);

    Partition<Count>* solidCounts = & dskGroup.getPartition<Count> ("solid");

    /** We configure the variant. */
    data.setSolid (solidCounts);

    /** We check that we got solid kmers. */
    if (solidCounts->getNbItems() == 0)  {  throw system::Exception ("This dataset has no solid kmers"); }

    DEBUG ((cout << "build_visitor : SortingCountAlgorithm END\n"));

    /** We save the state and kmer size at storage root level. */
    graph.getGroup().setProperty ("state",     Stringify::format("%d", graph._state));
    graph.getGroup().setProperty ("kmer_size", Stringify::format("%d", graph._kmerSize));
}


/* now build the rest of the graph */
template<typename Node, typename Edge, typename GraphDataVariant>
template <size_t span>  
void build_visitor_postsolid<Node,Edge,GraphDataVariant>::operator() (GraphData<span>& data) const 
{
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type  Type;

    /** We may have to stop just after configuration. */
    if (props->get(STR_CONFIG_ONLY))  { return; }

    if (!graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE))
    {
        throw system::Exception ("Graph construction failure during build_visitor_postsolid, the input _gatb/ folder (or .h5 file) needs to contain at least solid kmers");
    }

    size_t   kmerSize = graph.getKmerSize();

    // todo: see remark in build_visitor_solid, we should be able to get that info from elsewhere
    size_t minimizerSize = props->get(STR_MINIMIZER_SIZE)     ? props->getInt(STR_MINIMIZER_SIZE)      : 8;

    /** We create the kmer model. */
    data.setModel (new typename Kmer<span>::ModelCanonical (kmerSize));

    /** We may need a specific storage instance if we want to save the solid kmers in a separate file. */
    Storage* solidStorage = 0;
    if (props->get(STR_URI_SOLID_KMERS) != 0)
    {
        string solidsName = props->getStr(STR_URI_SOLID_KMERS);
        /** Historically (by convention), the file was deleted if its name is "none".
         * Now, since debloom may use the minimizer repartition function (stored in the solid file),
         * we must not delete the solid file. */
        bool autoDelete = false; // (solidsName == "none") || (solidsName == "null");
        solidStorage = StorageFactory(graph._storageMode).create (solidsName, false, autoDelete);
    }
    else
    {
        solidStorage = graph._storage;
    }
    LOCAL (solidStorage);

    /************************************************************/
    /*                         MPHF                             */
    /* note: theoretically could be done in parallel to debloom, but both tasks may or may not be IO intensive */
    /************************************************************/

    Group& dskGroup = (*solidStorage)("dsk"); 
    Partition<Count>* solidCounts = & dskGroup.getPartition<Count> ("solid");

    /** We create an instance of the MPHF Algorithm class (I was wondering: why is that a class, and not a function?) and execute it. */
    bool  noMphf = props->get("-no-mphf") != 0;
    if ((!noMphf) && (!graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE)))
    {
        DEBUG ((cout << "build_visitor : MPHFAlgorithm BEGIN\n"));

        /** We get the iterable for the solid counts and solid kmers. */
        Iterable<Type>*   solidKmers  = new IterableAdaptor<Count,Type,Count2TypeAdaptor<span> > (*solidCounts);

        MPHFAlgorithm<span> mphf_algo (
                dskGroup,
                "mphf",
                solidCounts,
                solidKmers,
                props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0, 
                // TODO enhancement: also pass the MAX_MEMORY parameter to enable or disable fast mode depending on it
                true  /* build=true, load=false */

                );
        graph.executeAlgorithm (mphf_algo, & graph.getStorage(), props, graph._info);
        data.setAbundance(mphf_algo.getAbundanceMap());
        data.setNodeState(mphf_algo.getNodeStateMap());
        data.setAdjacency(mphf_algo.getAdjacencyMap());
        graph.setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE);

        DEBUG ((cout << "build_visitor : MPHFAlgorithm END\n"));
    }

    /************************************************************/
    /*                         Bloom                            */
    /************************************************************/
    if (graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE) && !(graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BLOOM_DONE)))
    {
        DEBUG ((cout << "build_visitor : BloomAlgorithm BEGIN\n"));

        if (graph._bloomKind != BLOOM_NONE)
        {
            BloomAlgorithm<span> bloomAlgo (
                    graph.getStorage(),
                    data._solid,
                    kmerSize,
                    DebloomAlgorithm<span>::getNbBitsPerKmer (kmerSize, graph._debloomKind),
                    props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0,
                    graph._bloomKind
                    );
            graph.executeAlgorithm (bloomAlgo, & graph.getStorage(), props, graph._info);
            graph.setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BLOOM_DONE);
        }

        DEBUG ((cout << "build_visitor : BloomAlgorithm END\n"));
    }

    /************************************************************/
    /*                         Debloom                          */
    /************************************************************/
    if (graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BLOOM_DONE) && !(graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_DEBLOOM_DONE)))
    {
        DEBUG ((cout << "build_visitor : DebloomAlgorithm BEGIN\n"));

        Group& minimizersGroup = (graph.getStorage())("minimizers");

        /** We create a debloom instance and execute it. */
        DebloomAlgorithm<span>* debloom = DebloomAlgorithmFactory<span>::create (
                graph._debloomImpl,
                (graph.getStorage())("bloom"),
                (graph.getStorage())("debloom"),
                data._solid,
                kmerSize, minimizerSize,
                props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 0,
                props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0,
                graph._bloomKind,
                graph._debloomKind,
                "", 0,
                &minimizersGroup
                );
        LOCAL (debloom);

        graph.executeAlgorithm (*debloom, & graph.getStorage(), props, graph._info);

        graph.setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_DEBLOOM_DONE);

        /** We configure the variant. */
        data.setContainer (debloom->getContainerNode());

        DEBUG ((cout << "build_visitor : DebloomAlgorithm END\n"));
    }

    /************************************************************/
    /*                         Branching                        */
    /************************************************************/
    if (graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_DEBLOOM_DONE) && !(graph.checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BRANCHING_DONE)))
    {
        DEBUG ((cout << "build_visitor : BranchingAlgorithm BEGIN\n"));

        if (graph._branchingKind != BRANCHING_NONE)
        {
            BranchingAlgorithm<span, Node, Edge, GraphTemplate<Node,Edge,GraphDataVariant> > branchingAlgo (
                    graph,
                    graph.getStorage(),
                    graph._branchingKind,
                    props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0,
                    props
                    );
            graph.executeAlgorithm (branchingAlgo, & graph.getStorage(), props, graph._info);

            graph.setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_BRANCHING_DONE);

            /** We configure the variant. */
            data.setBranching (branchingAlgo.getBranchingCollection());
        }

        DEBUG ((cout << "build_visitor : BranchingAlgorithm END\n"));

    }


    /************************************************************/
    /*                    Post processing                       */
    /************************************************************/

    /** In case we choose another storage for the solid kmers, we have to update the graph state. */
    if (props->get(STR_URI_SOLID_KMERS) != 0)
    {
        data.setSolid (0);
        graph.unsetState (GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE);
    }

    /** We save library information in the root of the storage. */
    graph.getGroup().setProperty ("xml", string("\n") + LibraryInfo::getInfo().getXML());


    /** We save the state and kmer size at storage root level. */
    graph.getGroup().setProperty ("state",     Stringify::format("%d", graph._state));
    graph.getGroup().setProperty ("kmer_size", Stringify::format("%d", graph._kmerSize));

    /************************************************************/
    /*                        Clean up                          */
    /************************************************************/
}



/********************************************************************************
                 #####   ######      #     ######   #     #
                #     #  #     #    # #    #     #  #     #
                #        #     #   #   #   #     #  #     #
                #  ####  ######   #     #  ######   #######
                #     #  #   #    #######  #        #     #
                #     #  #    #   #     #  #        #     #
                 #####   #     #  #     #  #        #     #
********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
IOptionsParser* GraphTemplate<Node, Edge, GraphDataVariant>::getOptionsParser (bool includeMandatory)
{

    /** We build the root options parser. */
    OptionsParser* parser = new OptionsParser ("graph");

    /** We add children parser to it (kmer count, bloom/debloom, branching). */
    parser->push_back (SortingCountAlgorithm<>::getOptionsParser(includeMandatory));
    parser->push_back (DebloomAlgorithm<>::getOptionsParser());
    parser->push_back (BranchingAlgorithm<>::getOptionsParser());
    parser->push_front (new OptionNoParam  ("-no-mphf",       "don't construct the MPHF"));

    /** We create a "general options" parser. */
    IOptionsParser* parserGeneral  = new OptionsParser ("general");
    parserGeneral->push_front (new OptionOneParam (STR_INTEGER_PRECISION, "integers precision (0 for optimized value)", false, "0", false));
    parserGeneral->push_front (new OptionOneParam (STR_VERBOSE,           "verbosity level",      false, "1"  ));
    parserGeneral->push_front (new OptionOneParam (STR_NB_CORES,          "number of cores",      false, "0"  ));
    parserGeneral->push_front (new OptionNoParam  (STR_CONFIG_ONLY,       "dump config only"));
    
    parser->push_back  (parserGeneral);

    OptionsParser* parserDebug = new OptionsParser ("debug ");

    // those are only valid for GraphUnitigs, but GraphUnitigs doesn't have custom options (yet) so i'm adding here
    parserDebug->push_front (new OptionOneParam ("-nb-glue-partitions",       "number of glue partitions (automatically calculated by default)", false, "0"));
    //parserDebug->push_front (new OptionNoParam  ("-rebuild-graph",       "rebuild the whole graph starting from counted kmers"));
    parserDebug->push_front (new OptionNoParam  ("-skip-links",       "same, but       skip     links"));
    parserDebug->push_front (new OptionNoParam  ("-redo-links",       "same, but       redo     links"));
    parserDebug->push_front (new OptionNoParam  ("-skip-bglue",       "same, but       skip     bglue"));
    parserDebug->push_front (new OptionNoParam  ("-redo-bglue",       "same, but       redo     bglue"));
    parserDebug->push_front (new OptionNoParam  ("-skip-bcalm",       "same, but       skip     bcalm"));
    parserDebug->push_front (new OptionNoParam  ("-redo-bcalm",       "debug function, redo the bcalm algo"));

    /** We add it to the root parser. */
    parser->push_back  (parserDebug);

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
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>  GraphTemplate<Node, Edge, GraphDataVariant>::create (bank::IBank* bank, const char* fmt, ...)
{
    IOptionsParser* parser = getOptionsParser (false);   LOCAL(parser);

    /** We build the command line from the format and the ellipsis. */
    va_list args;
    va_start (args, fmt);
    std::string commandLine = Stringify::format(fmt, args);
    va_end (args);

    try
    {
        return  GraphTemplate(bank, parser->parseString(commandLine));
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
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>  GraphTemplate<Node, Edge, GraphDataVariant>::create (const char* fmt, ...)
{
    IOptionsParser* parser = getOptionsParser (true);   LOCAL (parser);

    /** We build the command line from the format and the ellipsis. */
    va_list args;
    va_start (args, fmt);
    std::string commandLine = Stringify::format(fmt, args);
    va_end (args);

    try
    {
        return  GraphTemplate (parser->parseString(commandLine)); /* will call the GraphTemplate<Node, Edge, GraphDataVariant>::GraphTemplate (tools::misc::IProperties* params) constructor */
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
template<typename Node, typename Edge, typename GraphDataVariant_t>
GraphTemplate<Node, Edge, GraphDataVariant_t>::GraphTemplate (size_t kmerSize)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant_t()), _kmerSize(kmerSize), _info("graph"),
      _state(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_INIT_DONE),
      _bloomKind(BLOOM_DEFAULT), _debloomKind(DEBLOOM_DEFAULT), _debloomImpl(DEBLOOM_IMPL_DEFAULT),
      _branchingKind(BRANCHING_STORED)
{
    /** We configure the data variant according to the provided kmer size. */
    setVariant (_variant, _kmerSize);

    /** We configure the graph data from the storage content. */
    boost::apply_visitor (configure_visitor<Node, Edge, GraphDataVariant_t>(*this, getStorage()),  *(GraphDataVariant_t*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE : loads a graph from a h5 file name or a _gatb/ folder
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant_t>
GraphTemplate<Node, Edge, GraphDataVariant_t>::GraphTemplate (const std::string& uri)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant_t()), _kmerSize(0), _info("graph"), 
      _name(System::file().getBaseName(uri))

{
    /** We create a storage instance. */
    /* (this is actually loading, not creating, the storage at "uri") */
    setStorage (StorageFactory(_storageMode).create (uri, false, false));

    /** We get some properties. */
    _state     = (GraphTemplate<Node, Edge, GraphDataVariant_t>::StateMask) atol (getGroup().getProperty ("state").c_str());
    _kmerSize  =                    atol (getGroup().getProperty ("kmer_size").c_str());

    /** We get library information in the root of the storage. */
    string xmlString = getGroup().getProperty ("xml");
    stringstream ss; ss << xmlString;   IProperties* props = new Properties(); LOCAL(props);
    props->readXML (ss);  getInfo().add (1, props);

    /** We configure the data variant according to the provided kmer size. */
    setVariant (_variant, _kmerSize);

    /** We configure the graph data from the storage content. */
    boost::apply_visitor (configure_visitor<Node, Edge, GraphDataVariant_t>(*this, getStorage()),  *(GraphDataVariant_t*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE : creates a graph from a bank
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>::GraphTemplate (bank::IBank* bank, tools::misc::IProperties* params)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph"),
      _state(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_INIT_DONE)
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    size_t integerPrecision = params->getInt (STR_INTEGER_PRECISION);

    /** We get other user parameters. */
    parse (params->getStr(STR_BLOOM_TYPE),        _bloomKind);
    parse (params->getStr(STR_DEBLOOM_TYPE),      _debloomKind);
    parse (params->getStr(STR_DEBLOOM_IMPL),      _debloomImpl);
    parse (params->getStr(STR_BRANCHING_TYPE),    _branchingKind);

    /** We configure the data variant according to the provided kmer size. */
    setVariant (_variant, _kmerSize, integerPrecision);

    /** We build the graph according to the wanted precision. */
    boost::apply_visitor (build_visitor_solid<Node, Edge, GraphDataVariant>(*this, bank,params),  *(GraphDataVariant*)_variant);
    boost::apply_visitor (build_visitor_postsolid<Node, Edge, GraphDataVariant>(*this, params),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE : creates (or completes; new feature) a graph from parsed command line arguments.
** INPUT   : a bank or a h5 file (new feature)
** OUTPUT  :
** RETURN  :
** REMARKS : this code could also work for (and is more generic than) the function above, TODO refactor.
            it is called from Graph::create (const char* fmt, ...)
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>::GraphTemplate (tools::misc::IProperties* params)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph"),
      _state(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_INIT_DONE)
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    size_t integerPrecision = params->getInt (STR_INTEGER_PRECISION);

    /** We get other user parameters. */
    parse (params->getStr(STR_BLOOM_TYPE),        _bloomKind);
    parse (params->getStr(STR_DEBLOOM_TYPE),      _debloomKind);
    parse (params->getStr(STR_DEBLOOM_IMPL),      _debloomImpl);
    parse (params->getStr(STR_BRANCHING_TYPE),    _branchingKind);

    /** We configure the data variant according to the provided kmer size. */
    setVariant (_variant, _kmerSize, integerPrecision);

    string input = params->getStr(STR_URI_INPUT);

    bool load_from_hdf5 = (system::impl::System::file().getExtension(input) == "h5");
    bool load_from_file = (system::impl::System::file().isFolderEndingWith(input,"_gatb"));
    bool load_graph = (load_from_hdf5 || load_from_file);
    if (load_graph)
    {
        /* it's not a bank, but rather a h5 file (kmercounted or more), let's complete it to a graph */
        
        cout << "Input is h5 or _gatb/ (we assume that kmer counting has already been done), we will complete it into a graph if necessary.\n"; 
        
        /** We create a storage instance. */
        /* (this is actually loading, not creating, the storage at "uri") */
        _storageMode = load_from_hdf5 ? STORAGE_HDF5 : STORAGE_FILE;
        bool append = true; // special storagehdf5 which will open the hdf5 file as read&write
        setStorage (StorageFactory(_storageMode).create (input, false, false, false, append));
    
        /** We get some properties. */
        _state     = (typename GraphTemplate<Node, Edge, GraphDataVariant>::StateMask) atol (getGroup().getProperty ("state").c_str());
        _kmerSize  =                    atol (getGroup().getProperty ("kmer_size").c_str());

        if (_kmerSize == 0) /* try the dsk group; this assumes kmer counting is done */
            _kmerSize  =    atol (getGroup("dsk").getProperty ("kmer_size").c_str());
        // also assume kmer counting is done
        setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_SORTING_COUNT_DONE);
    
        
        /* // doesn't work now with hdf5, because HDF5 attributes already exist and it doesn't like to overwrite them
        bool debug_rebuild= params->get("-rebuild-graph"); // this is a hidden debug option for now
        if (debug_rebuild)
            _state = 7; // init done, configure done, sorting_count done, but that's it
        */
        
        /** We get library information in the root of the storage. */
        string xmlString = getGroup().getProperty ("xml");
        stringstream ss; ss << xmlString;   IProperties* props = new Properties(); LOCAL(props);
        props->readXML (ss);  getInfo().add (1, props);
        
        /** We configure the data variant according to the provided kmer size. */
        setVariant (_variant, _kmerSize);

        /* call the configure visitor to load everything (e.g. solid kmers, MPHF, etc..) that's been done so far */
        boost::apply_visitor (configure_visitor<Node, Edge, GraphDataVariant>(*this, getStorage()),  *(GraphDataVariant*)_variant);

        boost::apply_visitor (build_visitor_postsolid<Node, Edge, GraphDataVariant>(*this, params),  *(GraphDataVariant*)_variant);
    }
    else
    {
        /** We build a Bank instance for the provided reads uri. */
        bank::IBank* bank = Bank::open (params->getStr(STR_URI_INPUT));

        /** We build the graph according to the wanted precision. */
        boost::apply_visitor (build_visitor_solid<Node, Edge, GraphDataVariant>(*this, bank,params),  *(GraphDataVariant*)_variant);
        boost::apply_visitor (build_visitor_postsolid<Node, Edge, GraphDataVariant>(*this, params),  *(GraphDataVariant*)_variant);
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
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>::GraphTemplate ()
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph"),
      _state(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_INIT_DONE),
      _bloomKind(BLOOM_DEFAULT),
      _debloomKind(DEBLOOM_DEFAULT), _debloomImpl(DEBLOOM_IMPL_DEFAULT), _branchingKind(BRANCHING_STORED)
{
    //std::cout << "empty graphtemplate constructor" << std::endl;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>::GraphTemplate (const GraphTemplate<Node, Edge, GraphDataVariant>& graph)
    : _storageMode(graph._storageMode), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(graph._kmerSize), _info("graph"), _name(graph._name), _state(graph._state)
{
    setStorage (graph._storage);

    if (graph._variant)  {  *((GraphDataVariant*)_variant) = *((GraphDataVariant*)graph._variant);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>& GraphTemplate<Node, Edge, GraphDataVariant>::operator= (const GraphTemplate<Node, Edge, GraphDataVariant>& graph)
{
    if (this != &graph)
    {
        _kmerSize        = graph._kmerSize;
        _storageMode     = graph._storageMode;
        _name            = graph._name;
        _info            = graph._info;
        _bloomKind       = graph._bloomKind;
        _debloomKind     = graph._debloomKind;
        _debloomImpl     = graph._debloomImpl;
        _branchingKind   = graph._branchingKind;
        _state           = graph._state;

        setStorage (graph._storage);

        if (graph._variant)  {  *((GraphDataVariant*)_variant) = *((GraphDataVariant*)graph._variant);  }
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
template<typename Node, typename Edge, typename GraphDataVariant>
GraphTemplate<Node, Edge, GraphDataVariant>::~GraphTemplate<Node, Edge, GraphDataVariant> ()
{
    //std::cout <<"normal graph destructor called" << std::endl;
    /** We release resources. */
    setStorage (0);
    if (_variant)  {  delete (GraphDataVariant*)_variant;  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::remove ()
{
    getStorage().remove();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : it's a bit of a misnomer, but in GATB language, Branching means "anything not simple". I.e. simple means outdegree=indegree=1 and branching is every other degree combination, including 0's
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
bool GraphTemplate<Node, Edge, GraphDataVariant>::isBranching (Node& node) const
{
    size_t in, out;
    degree(node, in, out); 
    return (! (in==1 && out==1));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
bool GraphTemplate<Node, Edge, GraphDataVariant>::isSimple (Edge& edge) const
{
    return this->outdegree(edge.from)==1  &&  this->indegree(edge.to)==1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
Node GraphTemplate<Node, Edge, GraphDataVariant>::reverse (const Node& node) const
{
    Node result = node;
    result.strand = kmer::StrandReverse(node.strand);
    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
BranchingNode_t<Node> GraphTemplate<Node, Edge, GraphDataVariant>::reverse (const BranchingNode_t<Node>& node) const
{
    BranchingNode_t<Node> result = node;
    result.strand = kmer::StrandReverse(node.strand);
    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
Edge GraphTemplate<Node, Edge, GraphDataVariant>::reverse (const Edge& edge) const
{
    Edge result;

    Nucleotide NT = edge.nt;

    if (edge.direction == DIR_INCOMING)
    {
        if (edge.from.strand == kmer::STRAND_FORWARD)  {  NT = (Nucleotide) edge.from.kmer[0];  }
        else                                           {  NT = kmer::reverse ((Nucleotide) edge.from.kmer[getKmerSize()-1]);  }
    }
    else
    {
        if (edge.from.strand == STRAND_FORWARD)  {  NT = kmer::reverse ((Nucleotide) edge.from.kmer[getKmerSize()-1]);  }
        else                                     {  NT = (Nucleotide) edge.from.kmer[0];                                }
    }

    result.set (
        edge.to.kmer,   edge.to.strand,
        edge.from.kmer, edge.from.strand,
        NT,
        edge.direction==DIR_OUTCOMING ? DIR_INCOMING : DIR_OUTCOMING
    );

    return result;
}



/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
size_t GraphTemplate<Node, Edge, GraphDataVariant>::indegree  (Node& node) const  {  return degree(node, DIR_INCOMING);   }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
size_t GraphTemplate<Node, Edge, GraphDataVariant>::outdegree (Node& node) const  {  return degree(node, DIR_OUTCOMING);  }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
size_t GraphTemplate<Node, Edge, GraphDataVariant>::degree (Node& node, Direction dir) const  {  return countNeighbors(node, dir);  } // used to be getNodes(node,dir).size() but made it faster

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::degree (Node& node, size_t &in, size_t &out) const  {  countNeighbors(node, in, out);  } 

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename Item, typename Functor, typename GraphDataVariant>
struct getItems_visitor : public boost::static_visitor<GraphVector<Item> >    {

    Node& source;  Direction direction;  Functor fct; 
    bool hasAdjacency;

    getItems_visitor (Node& aSource, Direction aDirection, bool hasAdjacency, Functor aFct) : source(aSource), direction(aDirection), fct(aFct), hasAdjacency(hasAdjacency) {}

    template<size_t span>  GraphVector<Item> operator() (const GraphData<span>& data) const
    {

        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        GraphVector<Item> items;
        GraphVector<Item> itemsAdj;

        size_t idx = 0;
        size_t idxAdj = 0;

        /** We get the specific typed value from the generic typed value. */
        Type sourceVal = source.template getKmer<Type>();

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();

        /* the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node */
        Type graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, kmerSize) );
        /* TODO opt: in some cases we may skip computing graine, e.g. whenever there are no neighbors */

        bool debug = false;
        /* use adjacency information when available, because it's faster than bloom */
        if (hasAdjacency)
        {
            unsigned long hashIndex = getNodeIndex<span>(data, source);
			if(hashIndex == ULLONG_MAX) return itemsAdj;

            unsigned char &value = (*(data._adjacency)).at(hashIndex);

            bool forwardStrand = (source.strand == STRAND_FORWARD);

            unsigned char bitmask;
         
 
            if (direction & DIR_OUTCOMING)
            {
                if (forwardStrand)
                    bitmask = value & 0xF;
                else
                {
                    bitmask = (value >> 4) & 0xF;

                    /* also revcomp the nt's: instead of GTCA (high bits to low), make it CAGT */
                    bitmask = ((bitmask & 3) << 2) | ((bitmask >> 2) & 3);
                }

                //std::cout << "getItems OUTCOMING: forward strand?" << forwardStrand << "; adjacency value: " << to_string((int)value) << " hashindex " << hashIndex << " dir " << direction << " bitmask " << (int)bitmask << std::endl;
           
                for (u_int64_t nt=0; nt<4; nt++)
                {
                    if (bitmask & (1 << nt))
                    {
                        typename Node::Value dest_value;
                        Strand dest_strand = STRAND_FORWARD;
                        Nucleotide dest_nt = (Nucleotide)nt;
                        Type forward, reverse;
                        Type single_nt; single_nt.setVal(nt);

                        forward = ( (graine << 2 )  + single_nt ) & mask; /* for speedup's, there is much we could precompute here */
                        reverse = revcomp (forward, kmerSize);

                        if (forward < reverse)
                            dest_value = forward;
                        else
                        {
                            dest_value = reverse;
                            dest_strand = STRAND_REVCOMP;
                        }

                        if (debug) std::cout << "adj, found OUT " << (dest_strand==STRAND_REVCOMP ? "REV" : "FWD") << " nt=" << nt << std::endl;
                        fct (itemsAdj, idxAdj++, source.kmer, source.strand, dest_value, dest_strand, dest_nt, DIR_OUTCOMING);
                    }
                }
            }

            if (direction & DIR_INCOMING)
            {
                if (forwardStrand)
                    bitmask = (value >> 4) & 0xF;
                else
                {
                    bitmask = value & 0xF;

                    /* also revcomp the nt's: instead of GTCA (high bits to low), make it CAGT */
                    bitmask = ((bitmask & 3) << 2) | ((bitmask >> 2) & 3);
                }

                //std::cout << "getItems INCOMING: forward strand?" << forwardStrand << "; adjacency value: " << to_string((int)value) << " hashindex " << hashIndex << " dir " << direction << " bitmask " << (int)bitmask << std::endl;

                /** IMPORTANT !!! Since we have hugely shift the nt value, we make sure to use a long enough integer. */
                for (u_int64_t nt=0; nt<4; nt++)
                {

                    if (bitmask & (1 << nt))
                    {
                        typename Node::Value dest_value;
                        Strand dest_strand = STRAND_FORWARD;
                        Nucleotide dest_nt = (Nucleotide)nt;
                        Type forward, reverse;
                        Type single_nt; single_nt.setVal(nt);

                        forward = ( (graine >> 2 )  + ( single_nt << ((kmerSize-1)*2)) ) & mask; /* previous kmer */
                        reverse = revcomp (forward, kmerSize);

                        if (forward < reverse)
                            dest_value = forward;
                        else
                        {
                            dest_value = reverse;
                            dest_strand = STRAND_REVCOMP;
                        }

                        if (debug) std::cout << "adj, found INC " << (dest_strand==STRAND_REVCOMP ? "REV" : "FWD") << " nt=" << nt << std::endl;
                        fct (itemsAdj, idxAdj++, source.kmer, source.strand, dest_value, dest_strand, dest_nt, DIR_INCOMING);
                    }
                }
            }

            /** We update the size of the container according to the number of found items. */
            itemsAdj.resize (idxAdj);

            /** We return the result. */
            return itemsAdj;
        }

        /* else, run classical neighbor queries using the data.contains() operation (bloom filters behind the scenes) */

        if (direction & DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                typename Node::Value dest_value;
                Type forward = ( (graine << 2 )  + nt) & mask;
                Type reverse = revcomp (forward, kmerSize);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        dest_value = forward;
                        if (debug) std::cout << "kmer  "<< sourceVal << " found OUT FWD nt=" << nt << std::endl;
                        fct (items, idx++, source.kmer, source.strand, dest_value, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        dest_value = reverse;
                        if (debug) std::cout << "found OUT REV nt=" << nt << std::endl;
                        fct (items, idx++, source.kmer, source.strand, dest_value, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
            }
        }

        if (direction & DIR_INCOMING)
        {
            /** IMPORTANT !!! Since we have hugely shift the nt value, we make sure to use a long enough integer. */
            for (u_int64_t nt=0; nt<4; nt++)
            {
                Type single_nt;
                single_nt.setVal(nt);
                single_nt <<=  ((kmerSize-1)*2);
                typename Node::Value dest_value;
                Type forward = ((graine >> 2 )  + single_nt ) & mask; /* previous kmer */
                Type reverse = revcomp (forward, kmerSize);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        dest_value = forward;
                        if (debug) std::cout << "found INC FWD nt=" << nt << std::endl;
                        // It used to be that "nt" was source[0] (in forward strand) and reverse(source[k-1]) in reverse, but i think it was wrong, so i changed it. TODO: delete this line if neighbors(..,DIR_INCOMING) causes no trouble for anyone, as it is now.
                        fct (items, idx++, source.kmer, source.strand, dest_value, STRAND_FORWARD, (Nucleotide)nt, DIR_INCOMING);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        dest_value = reverse;
                        if (debug) std::cout << "found INC REV nt=" << nt << std::endl;
                        fct (items, idx++, source.kmer, source.strand, dest_value, STRAND_REVCOMP, (Nucleotide)nt, DIR_INCOMING);
                    }
                }
            }
        }

        /** We update the size of the container according to the number of found items. */
        items.resize (idx);

        /** We return the result. */
        return items;
    }
};

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct Functor_getEdges {   void operator() (
    GraphVector<Edge>& items,
    size_t               idx,
    const typename Node::Value&   kmer_from,
    kmer::Strand         strand_from,
    const typename Node::Value&   kmer_to,
    kmer::Strand         strand_to,
    kmer::Nucleotide     nt,
    Direction            dir
) const
{
    items[idx++].set (kmer_from, strand_from, kmer_to, strand_to, nt, dir);
}};

/* TODO: so in principle, when Node is a NodeFast, we should be able to call the getItems_visitor 
 * operator directly, without apply_visitor (which seems to be expensive in the minia profiling using valgrind) 
 * but this seems tricky to code, so I gave up quickly */
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<Edge> GraphTemplate<Node, Edge, GraphDataVariant>::getEdges (Node source, Direction direction)  const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    return boost::apply_visitor (getItems_visitor<Node, Edge, Edge, Functor_getEdges<Node, Edge, GraphDataVariant>, GraphDataVariant>(source, direction, hasAdjacency, Functor_getEdges<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant);
}

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct Functor_getNodes {  void operator() (
    GraphVector<Node>&   items,
    size_t               idx,
    const typename Node::Value&   kmer_from,
    kmer::Strand         strand_from,
    const typename Node::Value&   kmer_to,
    kmer::Strand         strand_to,
    kmer::Nucleotide     nt,
    Direction            dir
) const
{
    items[idx++].set (kmer_to, strand_to);
}};

template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<Node> GraphTemplate<Node, Edge, GraphDataVariant>::getNodes (Node &source, Direction direction)  const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    return boost::apply_visitor (getItems_visitor<Node, Edge, Node, Functor_getNodes<Node, Edge, GraphDataVariant>, GraphDataVariant >(source, direction, hasAdjacency, Functor_getNodes<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : simple version of getITems above, just for counting number of neighbors. sorry for code duplication, I didn't want to use more functors (yet)
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct countNeighbors_visitor : public boost::static_visitor<void>    {

    Node& source;  Direction direction;
    bool hasAdjacency;
    size_t &indegree;
    size_t &outdegree;

    countNeighbors_visitor (Node& aSource, Direction aDirection, bool hasAdjacency, size_t& indegree, size_t& outdegree) : source(aSource), direction(aDirection), hasAdjacency(hasAdjacency),
    indegree(indegree), outdegree(outdegree) {}

    template<size_t span>  void operator() (const GraphData<span>& data) const
    {

        indegree = 0; outdegree = 0;

        //static int bitmask2nbneighbors[16] = { 
        
        /* use adjacency information when available, because it's faster than bloom */
        if (hasAdjacency)
        {
            unsigned long hashIndex = getNodeIndex<span>(data, source);
			if(hashIndex == ULLONG_MAX) return; // node was not found in the mphf 

            unsigned char &value = (*(data._adjacency)).at(hashIndex);

            bool forwardStrand = (source.strand == STRAND_FORWARD);

            unsigned char bitmask;
             
            if (forwardStrand)
                bitmask = value & 0xF;
            else
                bitmask = (value >> 4) & 0xF;

            outdegree = __builtin_popcount(bitmask);

            indegree = __builtin_popcount(value) - outdegree;

            return;
        }

        /* else, run classical neighbor queries using the data.contains() operation (bloom filters behind the scenes) */

        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        /** We get the specific typed value from the generic typed value. */
        Type sourceVal = source.template getKmer<Type>();

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();


        /* the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node */
        Type graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, kmerSize) );
        /* TODO opt: in some cases we may skip computing graine, e.g. whenever there are no neighbors */


        if (direction & DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                Type forward = ( (graine << 2 )  + nt) & mask;
                Type reverse = revcomp (forward, kmerSize);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        outdegree++;
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        outdegree++;
                    }
                }
            }
        }

        if (direction & DIR_INCOMING)
        {
            /** IMPORTANT !!! Since we have hugely shift the nt value, we make sure to use a long enough integer. */
            for (u_int64_t nt=0; nt<4; nt++)
            {
                Type single_nt;
                single_nt.setVal(nt);
                single_nt <<=  ((kmerSize-1)*2);
                Type forward = ((graine >> 2 )  + single_nt ) & mask; /* previous kmer */
                Type reverse = revcomp (forward, kmerSize);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        indegree++;
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        indegree++;
                    }
                }
            }
        }
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
unsigned char GraphTemplate<Node, Edge, GraphDataVariant>::countNeighbors (Node &source, Direction direction)  const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    size_t in, out;
    boost::apply_visitor (countNeighbors_visitor<Node, Edge, GraphDataVariant >(source, direction, hasAdjacency, in, out),  *(GraphDataVariant*)_variant);
    size_t res = 0;
    if (direction & DIR_INCOMING) res += in;
    if (direction & DIR_OUTCOMING) res += out;
    return res;
}

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::countNeighbors (Node &source, size_t &in, size_t &out)  const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    boost::apply_visitor (countNeighbors_visitor<Node, Edge, GraphDataVariant >(source, DIR_END, hasAdjacency, in, out),  *(GraphDataVariant*)_variant);
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : Why is'nt this code using getItems? ah, I guess, because there's a common functor called for both nodes neighbors at the same time. strange.
*********************************************************************/
template<typename Node, typename Edge, typename Item, typename Functor, typename GraphDataVariant>
struct getItemsCouple_visitor : public boost::static_visitor<GraphVector<pair<Item,Item> > >    {

    const Node& node1; const Node& node2;  Direction direction; Functor functor;

    getItemsCouple_visitor (const Node& node1, const Node& node2, Direction aDirection, Functor aFct)
        : node1(node1), node2(node2),  direction(aDirection), functor(aFct) {}

    template<size_t span>  GraphVector<pair<Item,Item> > operator() (const GraphData<span>& data) const
    {
        typedef typename Kmer<span>::Type  Type;

        size_t idx = 0;
        GraphVector < pair<Item,Item> > items;

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();

        /** We get the specific typed value from the generic typed value. */
        const Type& val1 = node1.template getKmer<Type>();
        const Type& val2 = node2.template getKmer<Type>();

        /* the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node */
        Type graine1 = ((node1.strand == STRAND_FORWARD) ?  val1 :  revcomp (val1, kmerSize) );
        Type graine2 = ((node2.strand == STRAND_FORWARD) ?  val2 :  revcomp (val2, kmerSize) );

        if (direction & DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                Type forward1 = ( (graine1 << 2 )  + nt) & mask;
                Type forward2 = ( (graine2 << 2 )  + nt) & mask;

                Type reverse1 = revcomp (forward1, kmerSize);
                Type reverse2 = revcomp (forward2, kmerSize);

                bool isForwardMin1 = forward1 < reverse1;
                bool isForwardMin2 = forward2 < reverse2;

                if (isForwardMin1==true && isForwardMin2==true)
                {
                    if (data.contains (forward1) && data.contains (forward2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        typename Node::Value dest_value1, dest_value2;
                        dest_value1 = forward1; dest_value2 = forward2;
                        functor (p,
                            node1.kmer, node1.strand, dest_value1, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, dest_value2, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
                else if (isForwardMin1==true && isForwardMin2==false)
                {
                    if (data.contains (forward1) && data.contains (reverse2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        typename Node::Value dest_value1, dest_value2;
                        dest_value1 = forward1; dest_value2 = reverse2;
                        functor (p,
                            node1.kmer, node1.strand, dest_value1, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, dest_value2, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
                else if (isForwardMin1==false && isForwardMin2==true)
                {
                    if (data.contains (reverse1) && data.contains (forward2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        typename Node::Value dest_value1, dest_value2;
                        dest_value1 = reverse1; dest_value2 = forward2;
                        functor (p,
                            node1.kmer, node1.strand, dest_value1, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, dest_value2, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
                else if (isForwardMin1==false && isForwardMin2==false)
                {
                    if (data.contains (reverse1) && data.contains (reverse2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        typename Node::Value dest_value1, dest_value2;
                        dest_value1 = reverse1; dest_value2 = reverse2;
                        functor (p,
                            node1.kmer, node1.strand, dest_value1, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, dest_value2, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
            }
        }

        if (direction & DIR_INCOMING)
        {
            throw system::ExceptionNotImplemented();
        }

        /** We update the size of the container according to the number of found items. */
        items.resize (idx);

        /** We return the result. */
        return items;
    }
};

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct Functor_getNodesCouple {  void operator() (
    pair<Node,Node>&   items,
    const typename Node::Value&   kmer_from1,
    kmer::Strand         strand_from1,
    const typename Node::Value&   kmer_to1,
    kmer::Strand         strand_to1,
    kmer::Nucleotide     nt1,
    Direction            dir1,
    const typename Node::Value&   kmer_from2,
    kmer::Strand         strand_from2,
    const typename Node::Value&   kmer_to2,
    kmer::Strand         strand_to2,
    kmer::Nucleotide     nt2,
    Direction            dir2
) const
{
    items.first.set  (kmer_to1, strand_to1);
    items.second.set (kmer_to2, strand_to2);
}};

template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<std::pair<Node,Node> > GraphTemplate<Node, Edge, GraphDataVariant>::getNodesCouple (const Node& node1, const Node& node2, Direction direction) const
{
    return boost::apply_visitor (getItemsCouple_visitor<Node, Edge, Node, Functor_getNodesCouple<Node, Edge, GraphDataVariant>, GraphDataVariant>(node1, node2, direction, Functor_getNodesCouple<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant);
}

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct Functor_getEdgesCouple {  void operator() (
        pair<Edge,Edge>&     items,
        const typename Node::Value&   kmer_from1,
        kmer::Strand         strand_from1,
        const typename Node::Value&   kmer_to1,
        kmer::Strand         strand_to1,
        kmer::Nucleotide     nt1,
        Direction            dir1,
        const typename Node::Value&   kmer_from2,
        kmer::Strand         strand_from2,
        const typename Node::Value&   kmer_to2,
        kmer::Strand         strand_to2,
        kmer::Nucleotide     nt2,
        Direction            dir2
) const
{
    items.first.set  (kmer_from1, strand_from1, kmer_to1, strand_to1, nt1, dir1);
    items.second.set (kmer_from2, strand_from2, kmer_to2, strand_to2, nt2, dir2);
}};

template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<std::pair<Edge,Edge> > GraphTemplate<Node, Edge, GraphDataVariant>::getEdgesCouple (const Node& node1, const Node& node2, Direction direction) const
{
    return boost::apply_visitor (getItemsCouple_visitor<Node, Edge, Edge, Functor_getEdgesCouple<Node, Edge, GraphDataVariant>, GraphDataVariant >(node1, node2, direction, Functor_getEdgesCouple<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct buildNode_visitor : public boost::static_visitor<Node>    {

    const tools::misc::Data& data;  size_t offset;

    buildNode_visitor (const tools::misc::Data& aData, size_t aOffset) : data(aData), offset(aOffset)  {}

    template<size_t span>  Node operator() (const GraphData<span>& graphData) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::ModelCanonical::Kmer Kmer;

        Kmer kmer = graphData._model->getKmer(data, offset);
        
        typename Node::Value dest_value;
        dest_value = kmer.value();
        return Node (dest_value, kmer.forward()==kmer.value() ? STRAND_FORWARD : STRAND_REVCOMP);
    }
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
Node GraphTemplate<Node, Edge, GraphDataVariant>::buildNode (const tools::misc::Data& data, size_t offset)  const
{
    return boost::apply_visitor (buildNode_visitor<Node, Edge, GraphDataVariant>(data,offset),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
Node GraphTemplate<Node, Edge, GraphDataVariant>::buildNode (const char* sequence)  const
{
    Data data ((char*)sequence);

    return boost::apply_visitor (buildNode_visitor<Node, Edge, GraphDataVariant>(data,0),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<BranchingNode_t<Node> > GraphTemplate<Node, Edge, GraphDataVariant>::getBranchingNodeNeighbors (Node& source, Direction direction) const
{
    GraphVector<BranchingNode_t<Node> >  result;

    /** We get the neighbors of the source node. */
    GraphVector<Edge> neighbors = this->neighborsEdge (source, direction);

    /** We resize the result vector. */
    result.resize (neighbors.size());

    /** We loop over all the neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** We get a simple path iterator from the current neighbor. */
        GraphIterator<Edge> path = this->simplePathEdge (neighbors[i].to, direction);

        /** We iterate this simple path from the current neighbor. */
        for (path.first(); !path.isDone(); path.next())  {}

        /** Note the trick here: we get the current path node, even if the path iteration is done. */
        Node& last = path.item().to;

        /** We set the ith branching neighbor node. */
        result[i].set (last.kmer, last.strand);
    }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<BranchingEdge_t<Node,Edge> > GraphTemplate<Node, Edge, GraphDataVariant>::getBranchingEdgeNeighbors (Node& source, Direction direction) const
{
    GraphVector<BranchingEdge_t<Node,Edge> >  result;

    /** We get the neighbors of the source node. */
    GraphVector<Edge> neighbors = this->neighborsEdge (source, direction);

    /** We resize the result vector. */
    result.resize (neighbors.size());

    /** We loop over all the neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        DEBUG ((cout << "neighbor[" << i << "] " << this->toString(neighbors[i]) << endl));

        /** We get a simple path iterator from the current neighbor. */
        GraphIterator<Edge> path = this->simplePathEdge (neighbors[i].to, direction);

        /** We iterate this simple path from the current neighbor. */
        for (path.first(); !path.isDone(); path.next())
        {
            DEBUG ((cout << "===> " << this->toString(*path) << endl));
        }

        /** Note the trick here: we get the current path node, even if the path iteration is done. */
        Node& last = path.item().to;

        /** We set the ith branching neighbor node. By convention, we memorize the nucleotide of the initial transition. */
        result[i].set (source.kmer, source.strand, last.kmer, last.strand, neighbors[i].nt, direction, path.rank()+1);
    }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename Item, typename Functor>
struct getItem_visitor : public boost::static_visitor<Item>    {

    Node& source;  Direction direction;  Nucleotide nt; bool& exists; Functor fct; bool hasAdjacency;

    getItem_visitor (Node& aSource, Direction aDirection, Nucleotide aNt, bool hasAdjacency, bool& aExists, Functor aFct)
        : source(aSource), direction(aDirection), nt(aNt), exists(aExists), fct(aFct), hasAdjacency(hasAdjacency) {}

    template<size_t span>  Item operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        Item item;
        Item itemAdj;

        /** We get the specific typed value from the generic typed value. */
        const Type& sourceVal = source.template getKmer<Type>();

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();

        /* the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node */
        Type graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, kmerSize) );


        // this is basically the code in getItems_visitor without a loop on nt
        bool debug = false;
        if (hasAdjacency)
        {
            unsigned long hashIndex = getNodeIndex<span>(data, source);
			if(hashIndex == ULLONG_MAX) {exists = false; return itemAdj;} // node was not found in the mphf 

            unsigned char &value = (*(data._adjacency)).at(hashIndex);

            bool forwardStrand = (source.strand == STRAND_FORWARD);

            unsigned char bitmask;

            if (direction & DIR_OUTCOMING)
            {
                if (forwardStrand)
                    bitmask = value & 0xF;
                else
                {
                    bitmask = (value >> 4) & 0xF;

                    /* also revcomp the nt's: instead of GTCA (high bits to low), make it CAGT */
                    bitmask = ((bitmask & 3) << 2) | ((bitmask >> 2) & 3);
                }

                //std::cout << "getItems OUTCOMING: forward strand?" << forwardStrand << "; adjacency value: " << to_string((int)value) << " hashindex " << hashIndex << " dir " << direction << " bitmask " << (int)bitmask << std::endl;

                if (bitmask & (1 << nt))
                {
                    typename Node::Value dest_value;
                    Strand dest_strand = STRAND_FORWARD;
                    Nucleotide dest_nt = (Nucleotide)nt;
                    Type forward, reverse;
                    Type single_nt; single_nt.setVal(nt);

                    forward = ( (graine << 2 )  + single_nt ) & mask; /* for speedup's, there is much we could precompute here */
                    reverse = revcomp (forward, kmerSize);

                    if (forward < reverse)
                        dest_value = forward;
                    else
                    {
                        dest_value = reverse;
                        dest_strand = STRAND_REVCOMP;
                    }

                    if (debug) std::cout << "getItem adj, found OUT " << (dest_strand==STRAND_REVCOMP ? "REV" : "FWD") << " nt=" << nt << std::endl;
                    fct (item, source.kmer, source.strand, dest_value, dest_strand, dest_nt, DIR_OUTCOMING);
                    exists = true;
                }
                else
                    exists = false;
            }

            if (direction & DIR_INCOMING)
            {
                if (forwardStrand)
                    bitmask = (value >> 4) & 0xF;
                else
                {
                    bitmask = value & 0xF;

                    /* also revcomp the nt's: instead of GTCA (high bits to low), make it CAGT */
                    bitmask = ((bitmask & 3) << 2) | ((bitmask >> 2) & 3);
                }

                //std::cout << "getItems INCOMING: forward strand?" << forwardStrand << "; adjacency value: " << to_string((int)value) << " hashindex " << hashIndex << " dir " << direction << " bitmask " << (int)bitmask << std::endl;

                if (bitmask & (1 << nt))
                {
                    typename Node::Value dest_value;
                    Strand dest_strand = STRAND_FORWARD;
                    Nucleotide dest_nt = (Nucleotide)nt;
                    Type forward, reverse;
                    Type single_nt; single_nt.setVal(nt);

                    forward = ( (graine >> 2 )  + ( single_nt << ((kmerSize-1)*2)) ) & mask; /* previous kmer */
                    reverse = revcomp (forward, kmerSize);

                    if (forward < reverse)
                        dest_value = forward;
                    else
                    {
                        dest_value = reverse;
                        dest_strand = STRAND_REVCOMP;
                    }

                    if (debug) std::cout << "getItem adj, found INC " << (dest_strand==STRAND_REVCOMP ? "REV" : "FWD") << " nt=" << nt << std::endl;
                    fct (item, source.kmer, source.strand, dest_value, dest_strand, dest_nt, DIR_INCOMING);
                    exists = true;
                }
                else
                    exists = false;
            }

            /** We return the result. */
            return item;
        }

        if (direction & DIR_OUTCOMING)
        {
            typename Node::Value dest_value;
            Type forward = ( (graine << 2 )  + nt) & mask;
            Type reverse = revcomp (forward, kmerSize);

            if (forward < reverse)
            {
                if (data.contains (forward))
                {
                    dest_value = forward;
                    fct (item, source.kmer, source.strand, dest_value, STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                    exists = true;
                }
                else
                {
                    exists = false;
                }
            }
            else
            {
                if (data.contains (reverse))
                {
                    dest_value = reverse;
                    fct (item, source.kmer, source.strand, dest_value, STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
                    exists = true;
                }
                else
                {
                    exists = false;
                }
            }
        }

        if (direction & DIR_INCOMING)
        {
            Type single_nt;
            single_nt.setVal(nt);
            single_nt <<= ((kmerSize-1)*2);
            typename Node::Value dest_value;
            Type forward = ((graine >> 2 )  + single_nt ) & mask; /* previous kmer */
            Type reverse = revcomp (forward, kmerSize);

            if (forward < reverse)
            {
                if (data.contains (forward))
                {
                    dest_value = forward;
                    fct (item, source.kmer, source.strand, dest_value, STRAND_FORWARD, (Nucleotide)nt, DIR_INCOMING);
                    exists = true;
                }
                else
                {
                    exists = false;
                }
            }
            else
            {
                if (data.contains (reverse))
                {
                    dest_value = reverse;
                    fct (item, source.kmer, source.strand, dest_value, STRAND_REVCOMP, (Nucleotide)nt, DIR_INCOMING);
                    exists = true;
                }
                else
                {
                    exists = false;
                }
            }
        }

        /** We return the result. */
        return item;
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
struct Functor_getNode {  void operator() (
    Node&                item,
    const typename Node::Value&   kmer_from,
    kmer::Strand         strand_from,
    const typename Node::Value&   kmer_to,
    kmer::Strand         strand_to,
    kmer::Nucleotide     nt,
    Direction            dir
) const
{
    item.set (kmer_to, strand_to);
}};


template<typename Node, typename Edge, typename GraphDataVariant>
Node GraphTemplate<Node, Edge, GraphDataVariant>::getNode (Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    return boost::apply_visitor (getItem_visitor<Node, Edge, Node, Functor_getNode<Node, Edge, GraphDataVariant> >(source, dir, nt, hasAdjacency, exists, Functor_getNode<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<Edge> GraphTemplate<Node, Edge, GraphDataVariant>::getEdgeValues (const typename Node::Value& kmer) const
{
    Node source (kmer);

    GraphVector<Edge> v1 = getEdges(source,          DIR_OUTCOMING);
    GraphVector<Edge> v2 = getEdges(reverse(source), DIR_OUTCOMING);
#if 0
    v1.insert (v1.end(), v2.begin(), v2.end());
#else
    size_t n1=v1.size(), n2=v2.size();
    v1.resize (n1+n2);
    for (size_t i=0; i<n2; i++)  { v1[i+n1] = v2[i];  }
#endif

    return v1;
}

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<Node> GraphTemplate<Node, Edge, GraphDataVariant>::getNodeValues (const typename Node::Value& kmer) const
{
    Node source (kmer);
    Node rev_source (reverse(source));

    GraphVector<Node> v1 = getNodes(source,          DIR_OUTCOMING);
    GraphVector<Node> v2 = getNodes(rev_source, DIR_OUTCOMING);

#if 0
    v1.insert (v1.end(), v2.begin(), v2.end());
#else
    size_t n1=v1.size(), n2=v2.size();
    v1.resize (n1+n2);
    for (size_t i=0; i<n2; i++)  { v1[i+n1] = v2[i];  }
#endif

    return v1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<BranchingEdge_t<Node,Edge> > GraphTemplate<Node, Edge, GraphDataVariant>::getBranchingEdgeValues (const typename Node::Value& kmer) const
{
    Node source (kmer);
    Node rev_source (reverse(source));

    GraphVector<BranchingEdge_t<Node,Edge> > v1 = getBranchingEdgeNeighbors (source,          DIR_OUTCOMING);
    GraphVector<BranchingEdge_t<Node,Edge> > v2 = getBranchingEdgeNeighbors (rev_source, DIR_OUTCOMING);
#if 0
    v1.insert (v1.end(), v2.begin(), v2.end());
#else
    size_t n1=v1.size(), n2=v2.size();
    v1.resize (n1+n2);
    for (size_t i=0; i<n2; i++)  { v1[i+n1] = v2[i];  }
#endif

    return v1;
}

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphVector<BranchingNode_t<Node> > GraphTemplate<Node, Edge, GraphDataVariant>::getBranchingNodeValues (const typename Node::Value& kmer) const
{
    Node source (kmer);
    Node rev_source (reverse(source));

    GraphVector<BranchingNode_t<Node> > v1 = getBranchingNodeNeighbors (source,          DIR_OUTCOMING);
    GraphVector<BranchingNode_t<Node> > v2 = getBranchingNodeNeighbors (rev_source, DIR_OUTCOMING);


#if 0
    v1.insert (v1.end(), v2.begin(), v2.end());
#else
    size_t n1=v1.size(), n2=v2.size();
    v1.resize (n1+n2);
    for (size_t i=0; i<n2; i++)  { v1[i+n1] = v2[i];  }
#endif

    return v1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct contains_visitor : public boost::static_visitor<bool>    {

    const Node& node;
    contains_visitor (const Node& aNode) : node(aNode) {}

    template<size_t span>  bool operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        return data.contains (node.template getKmer<Type>());
    }
};

/********************************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
bool GraphTemplate<Node, Edge, GraphDataVariant>::contains (const Node& item) const
{
    return boost::apply_visitor (contains_visitor<Node, Edge, GraphDataVariant>(item),  *(GraphDataVariant*)_variant);
}

template<typename Node, typename Edge, typename GraphDataVariant>
template<size_t span>
bool GraphTemplate<Node, Edge, GraphDataVariant>::contains (const typename Kmer<span>::Type& item) const
{
    return boost::apply_visitor (contains_visitor<Node, Edge, GraphDataVariant>(item),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename NodeType, typename GraphDataVariant>  struct BranchingFilter
{
    const GraphTemplate<Node, Edge, GraphDataVariant>& graph;
    BranchingFilter(const GraphTemplate<Node, Edge, GraphDataVariant>& graph) : graph(graph) {}
    bool operator () (NodeType& item) { return graph.isBranching(item); }
};

template<typename Node, typename Edge, typename NodeType, typename GraphDataVariant>
struct nodes_visitor : public boost::static_visitor<tools::dp::ISmartIterator<NodeType>*>
{
    const GraphTemplate<Node, Edge, GraphDataVariant>& graph;
    nodes_visitor (const GraphTemplate<Node, Edge, GraphDataVariant>& graph) : graph(graph) {}

    template<size_t span>  tools::dp::ISmartIterator<NodeType>* operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::Count Count;

        // soo.. we're defining a iterator class inside a visitor. that's just what it is.
        class NodeIterator : public tools::dp::ISmartIterator<NodeType>
        {
        public:
            NodeIterator (tools::dp::Iterator<Count>* ref, u_int64_t nbItems)
                : _ref(0),  _rank(0), _isDone(true), _nbItems(nbItems)   {  
                    setRef(ref);  
                    this->_item->strand = STRAND_FORWARD;  // iterated nodes are always in forward strand.
                }

            ~NodeIterator ()  { setRef(0);   }

            u_int64_t rank () const { return _rank; }

            /** \copydoc  Iterator::first */
            void first()
            {
                _ref->first();
                _rank   = 0;
                _isDone = _ref->isDone();

                if (!_isDone)
                {
                    // NOTE: doesn't check if node is deleted (as it would be expensive to compute MPHF index)
                    this->_rank ++;
                    this->_item->kmer      = _ref->item().value;
                    this->_item->abundance = _ref->item().abundance;
                    this->_item->mphfIndex = 0;
                    this->_item->iterationRank = this->_rank;
                }
            }

            /** \copydoc  Iterator::next */
            void next()
            {
                _ref->next();
                _isDone = _ref->isDone();
                if (!_isDone)
                {
                    // NOTE: doesn't check if node is deleted (as it would be expensive to compute MPHF index)
                    this->_rank ++;
                    this->_item->kmer      = _ref->item().value;
                    this->_item->abundance = _ref->item().abundance;
                    this->_item->mphfIndex = 0;
                    this->_item->iterationRank = this->_rank;
                }
            }

            /** \copydoc  Iterator::isDone */
            bool isDone() { return _isDone;  }

            /** \copydoc  Iterator::item */
            NodeType& item ()  {  return *(this->_item);  }

            /** */
            void setItem (NodeType& i)
            {
                /** We set the node item to be set for the current iterator. */
                this->_item = &i;
                this->_item->strand = STRAND_FORWARD;

                /** We set the kmer item to be set for the kmer iterator. */
                // _ref->setItem (i.kmer.value.get<T>());
                // TODO doc: uh, why is that commented?
            }

            /** */
            u_int64_t size () const { return _nbItems; }

        private:
            tools::dp::Iterator<Count>* _ref;
            void setRef (tools::dp::Iterator<Count>* ref)  { SP_SETATTR(ref); }

            u_int64_t _rank;
            bool      _isDone;
            u_int64_t _nbItems;
        };

        // now this is the actual code for returning a node iterator, apparently

        if (typeid(NodeType) == typeid(Node) )
        {
            if (data._solid != 0)
            {
                return new NodeIterator (data._solid->iterator (), data._solid->getNbItems());
            }
            else
            {
                throw system::Exception("Iteration impossible (no solid nodes available)");
            }
        }
        else if (typeid(NodeType) == typeid(BranchingNode_t<Node>))
        {
            if (data._branching != 0)
            {
                /** We have a branching container*/
                return new NodeIterator (data._branching->iterator (), data._branching->getNbItems());
            }
            else if (data._solid != 0)
            {
                /** We don't have pre-computed branching nodes container. We have to compute them on the fly
                 * from the solid kmers. We can do that by filtering out all non branching nodes. */
                return new FilterIterator<NodeType,BranchingFilter<Node, Edge, NodeType, GraphDataVariant> > (
                    new NodeIterator (data._solid->iterator (), data._solid->getNbItems()),
                    BranchingFilter<Node, Edge, NodeType, GraphDataVariant> (graph)
                );
            }
            else
            {
                throw system::Exception("Iteration impossible (no solid nor branching nodes available)");
            }
        }
        else {  throw system::Exception("Invalid type");  }
    }
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphIterator<Node> GraphTemplate<Node, Edge, GraphDataVariant>::getNodes () const
{
    return GraphIterator<Node> (boost::apply_visitor (nodes_visitor<Node,Edge, Node, GraphDataVariant>(*this),  *(GraphDataVariant*)_variant));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
GraphIterator<BranchingNode_t<Node> > GraphTemplate<Node, Edge, GraphDataVariant>::getBranchingNodes () const
{
    return GraphIterator<BranchingNode_t<Node> > (boost::apply_visitor (nodes_visitor<Node, Edge, BranchingNode_t<Node>, GraphDataVariant>(*this),  *(GraphDataVariant*)_variant));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct toString_node_visitor : public boost::static_visitor<std::string>    {

    const Node& node;
    toString_node_visitor (const Node& aNode) : node(aNode) {}

    template<size_t span>  std::string operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        Type value = node.template getKmer<Type>();
        if (node.strand == STRAND_FORWARD)   {  return data._model->toString (value);  }
        else                                 {  return data._model->toString (data._model->reverse (value));  }
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::toString (const Node& node) const
{
    return boost::apply_visitor (toString_node_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct debugString_node_visitor : public boost::static_visitor<std::string>    {

    const Node& node;  kmer::Strand strand; int mode;
    debugString_node_visitor (const Node& aNode, kmer::Strand aStrand, int aMode) : node(aNode), strand(aStrand), mode(aMode) {}

    template<size_t span>  std::string operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        std::stringstream ss;

        /** We set the strings for the kmer and the strand. */
        std::string kmerStr;
        std::string strandStr;

        Type value = node.template getKmer<Type>();

        if (strand == STRAND_ALL || node.strand == strand)
        {
            kmerStr = data._model->toString (value);
            strandStr = (node.strand==STRAND_FORWARD ? "FWD" : "REV");
        }
        else
        {
            Type reverse = data._model->reverse (value);
            kmerStr = data._model->toString (reverse);
            strandStr = (node.strand==STRAND_FORWARD ? "REV" : "FWD");
        }

             if (mode==0)  {  ss << "[ " << value <<  " / " << strandStr << "]";  }
        else if (mode==1)  {  ss << "[ " << kmerStr <<  " / " << strandStr << "]";  }
        else if (mode==2)  {  ss << "[ " << kmerStr <<  " / " << strandStr << "]";  }

        return ss.str();

    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::debugString (const Node& node, kmer::Strand strand, int mode) const
{
    return boost::apply_visitor (debugString_node_visitor<Node, Edge, GraphDataVariant>(node,strand,mode),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
struct debugString_edge_visitor : public boost::static_visitor<std::string>    {

    const Edge& edge;  Strand strand;  int mode;
    debugString_edge_visitor (const Edge& aEdge, Strand aStrand, int aMode=0) : edge(aEdge), strand(aStrand), mode(aMode) {}

    template<size_t span>  std::string operator() (const GraphData<span>& data) const
    {
        std::stringstream ss;

        if (mode == 0)
        {
            ss << "["
               << "(" << edge.from.kmer << "," << edge.from.strand-1 << ")"
               << " ";
            if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
            else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }

            ss  << " "
               << "(" << edge.to.kmer << "," << edge.to.strand-1 << ")"
               << "]";
        }
        else if (mode==1)
        {
            ss << "["
               << debugString_node_visitor<Node, Edge, GraphDataVariant>(edge.from, strand, 1) (data)
               << " ";
            if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
            else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }

            ss  << " "
                    << debugString_node_visitor<Node, Edge, GraphDataVariant>(edge.to, strand, 1) (data)
               << "]";
        }
        else if (mode==2)
        {
            ss << "["
               << toString_node_visitor<Node, Edge, GraphDataVariant>(edge.from) (data)
               << " ";
            if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
            else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }

            ss  << " "
                    << toString_node_visitor<Node, Edge, GraphDataVariant>(edge.to) (data)
               << "]";
        }

        return ss.str();
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::debugString (const Edge& edge, kmer::Strand strand, int mode) const
{
    return boost::apply_visitor (debugString_edge_visitor<Node, Edge, GraphDataVariant>(edge, strand, mode),  *(GraphDataVariant*)_variant);
}

template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::toString (const Edge& edge) const
{
    std::stringstream ss;

    ss << "["  << this->toString (edge.from)  << " ";
    if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
    else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }
    ss << " "  << this->toString (edge.to)  << "]";

    return ss.str();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::toString (const BranchingEdge_t<Node,Edge>& edge) const
{
    std::stringstream ss;

    ss << "["  << this->toString (edge.from)  << " ";
    if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "," << edge.distance << "-->";  }
    else                                  {  ss <<  "<--" << ascii(edge.nt) << "," << edge.distance << "--";   }
    ss << " "  << this->toString (edge.to)  << "]";

    return ss.str();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
int GraphTemplate<Node, Edge, GraphDataVariant>::simplePathAvance (Node& node, Direction dir, kmer::Nucleotide& nt) const
{
    Edge edge;
    int res = simplePathAvance (node, dir, edge);
    nt = edge.nt;
    return res;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
int GraphTemplate<Node, Edge, GraphDataVariant>::simplePathAvance (Node& node, Direction dir) const
{
    Edge output;  return simplePathAvance (node, dir, output);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant>
int GraphTemplate<Node, Edge, GraphDataVariant>::simplePathAvance (Node& node, Direction dir, Edge& output) const
{
    GraphVector<Edge> neighbors = this->neighborsEdge (node, dir);

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

template<typename Node, typename Edge, typename GraphDataVariant>
double GraphTemplate<Node, Edge, GraphDataVariant>::
simplePathMeanAbundance     (Node& node, Direction dir) const
{
    GraphIterator <Node> itNodes = simplePath (node, dir);

    unsigned long mean_abundance = 0;
    int length = 0;
    for (itNodes.first(); !itNodes.isDone(); itNodes.next())
    {
        unsigned int abundance = queryAbundance(*itNodes);
        mean_abundance += abundance;
        length++;
    }
    if (length == 0)
        return 0;
    double meanAbundance = (double)mean_abundance / ((double)length);
    double stdevAbundance = 0;

    // get std dev, for debug only
    bool debugstdev = false;
    if (debugstdev)
    {
        for (itNodes.first(); !itNodes.isDone(); itNodes.next())
        {
            unsigned int abundance = queryAbundance((*itNodes));
            stdevAbundance += pow(fabs(abundance-meanAbundance),2);
        }
        stdevAbundance = sqrt(stdevAbundance/(double)length);
    }

    return meanAbundance;
}

template<typename Node, typename Edge, typename GraphDataVariant>
Node  GraphTemplate<Node, Edge, GraphDataVariant>::
simplePathLastNode  (Node& node, Direction dir) const
{
    GraphIterator <Node> itNodes = simplePath (node, dir);
    itNodes.first();
    Node cur = *itNodes;
    for (; !itNodes.isDone(); itNodes.next())
        cur = *itNodes; // is there an easier way?
    return cur;
}

template<typename Node, typename Edge, typename GraphDataVariant>
unsigned int GraphTemplate<Node, Edge, GraphDataVariant>::
simplePathLength (Node& node, Direction dir) const
{
    GraphIterator <Node> itNodes = simplePath (node, dir);
    unsigned int length = 0;
    for (itNodes.first(); !itNodes.isDone(); itNodes.next())
        length++;
    return length;
}

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::
unitigDelete (Node& node, Direction dir, NodesDeleter<Node,Edge, GraphTemplate<Node, Edge, GraphDataVariant>> &nodesDeleter) 
{
    std::cout << "Graph::unitigDelete not implemented" << std::endl; exit(1);
}

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::
unitigDelete (Node& node) 
{
    std::cout << "Graph::simplePathDelete not implemented" << std::endl; exit(1);
}


template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::
simplePathDelete          (Node& node, Direction dir, NodesDeleter<Node,Edge, GraphTemplate<Node, Edge, GraphDataVariant>>& nodesDeleter) 
{
    GraphIterator <Node> itNodes = simplePath (node, dir);
    for (itNodes.first(); !itNodes.isDone(); itNodes.next())
        nodesDeleter.markToDelete(*itNodes);
    nodesDeleter.markToDelete(node); // don't forget the start node
}

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::
unitigMark            (Node& node)  // used to flag simple path as traversed, in minia
{
    std::cout << "Graph::simplePathMark not implemented. Only in GraphUnitigs it is." << std::endl; exit(1);
}

template<typename Node, typename Edge, typename GraphDataVariant>
bool GraphTemplate<Node, Edge, GraphDataVariant>::
unitigIsMarked        (Node& node)
{
    std::cout << "Graph::simplePathIsMarked not implemented" << std::endl; exit(1);
}





template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::
unitigSequence (Node& node, bool& isolatedLeft, bool& isolatedRight) const
{
    std::cout << "Graph::simplePathSequence not implemented (only in GraphUnitigs yet)" << std::endl;// would be easy to implement, just see Minia's legacy code for example
    return "";
}

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::
simplePathLongest_avance(const Node& node, Direction dir, int& seqLength, int& endDegree, bool dummy, float& coverage, std::string* seq, std::vector<Node> *unitigNodes) 
{
    std::cout << "Graph::simplePathLongest_avance not implemented (only in GraphUnitigs yet)" << std::endl;// would be easy to implement, just see Minia's legacy code for example
    return ;
}

template<typename Node, typename Edge, typename GraphDataVariant>
std::string GraphTemplate<Node, Edge, GraphDataVariant>::
simplePathBothDirections(const Node& node, bool& isolatedLeft, bool& isolatedRight, bool dummy, float& coverage) 
{
    std::cout << "Graph::simplePathLongest not implemented (only in GraphUnitigs yet)" << std::endl;// if wanted to implement, just see Minia's legacy code for example
    return "";
}


template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::
debugPrintAllUnitigs() const
{
    std::cout << "Graph::debugPritnAllUnitigs implemented (only in GraphUnitigs yet)" << std::endl;
}



/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename Item, typename Functor, typename GraphDataVariant>
class AbstractSimplePathIterator : public tools::dp::ISmartIterator<Item>
{
public:

    AbstractSimplePathIterator (const GraphTemplate<Node, Edge, GraphDataVariant>& graph, const Node& node, Direction dir, const Functor& update)
        : _graph(graph), _dir(dir), _rank(0), _isDone(true), _update(update) {}

    virtual ~AbstractSimplePathIterator() {}

    u_int64_t rank () const { return _rank; }

    /** \copydoc  Iterator::first */
    void first()
    {
        _rank   = 0;
        _isDone = false;
        next ();
        /** The rank must begin at 0; we just remove the 1 added in the first call to next. */
        _rank --;
    }

    /** \copydoc  Iterator::next */
    void next()
    {
        _rank ++;
        _update (_graph, *(this->_item), _dir, _isDone);
    }

    /** \copydoc  Iterator::isDone */
    bool isDone() {  return _isDone;  }

    /** \copydoc  Iterator::item */
    Item& item ()  {  return *(this->_item);  }

    /** */
    u_int64_t size () const { return 0; }

protected:
    const GraphTemplate<Node, Edge, GraphDataVariant>&       _graph;
    Direction          _dir;
    u_int64_t          _rank;
    bool               _isDone;
    const Functor&     _update;
};

/** */
template<typename Node, typename Edge, typename Functor, typename GraphDataVariant>
class NodeSimplePathIterator : public AbstractSimplePathIterator<Node,Edge, Node, Functor, GraphDataVariant>
{
public:
    NodeSimplePathIterator (const GraphTemplate<Node, Edge, GraphDataVariant>& graph, const Node& node, Direction dir, const Functor& update)
        : AbstractSimplePathIterator<Node,Edge, Node,Functor, GraphDataVariant> (graph, node, dir, update)  { *(this->_item) = node;  }
};

/** */
template<typename Node, typename Edge, typename Functor, typename GraphDataVariant>
class EdgeSimplePathIterator : public AbstractSimplePathIterator<Node,Edge, Edge, Functor, GraphDataVariant>
{
public:
    EdgeSimplePathIterator (const GraphTemplate<Node, Edge, GraphDataVariant>& graph, const Node& node, Direction dir, const Functor& udpate)
        : AbstractSimplePathIterator<Node,Edge, Edge, Functor, GraphDataVariant> (graph, node, dir, udpate)
    { this->_item->set (node.kmer, node.strand, node.kmer, node.strand, NUCL_UNKNOWN, dir); }
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant> 
struct Functor_getSimpleNodeIterator {  void operator() (
    const GraphTemplate<Node, Edge, GraphDataVariant>&         graph,
    Node&                item,
    Direction            dir,
    bool&                isDone
) const
{
    Edge output;

    /** We check if we have a simple node. */
    if (graph.simplePathAvance (item, dir, output) > 0)
    {
        /** We update the currently iterated node. */
        item = output.to;
    }
    else
    {
        /** We can't have a non branching node in the wanted direction => iteration is finished. */
        isDone = true;
    }
}};

template<typename Node, typename Edge, typename GraphDataVariant> 
GraphIterator<Node> GraphTemplate<Node, Edge, GraphDataVariant>::getSimpleNodeIterator (Node& node, Direction dir) const
{
    return GraphIterator<Node> (new NodeSimplePathIterator <Node, Edge, Functor_getSimpleNodeIterator<Node, Edge, GraphDataVariant>, GraphDataVariant> (*this, node, dir, Functor_getSimpleNodeIterator<Node, Edge, GraphDataVariant>()));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant> 
struct Functor_getSimpleEdgeIterator {  void operator() (
    const GraphTemplate<Node, Edge, GraphDataVariant>&         graph,
    Edge&                item,
    Direction            dir,
    bool&                isDone
) const
{
    Edge output;

    /** We check if we have a simple node. */
    int res = graph.simplePathAvance (item.to, dir, output);

    /** NOTE: we update the item in case we have outdegree==1 (case>0 and case==-2) */
    if (res > 0  ||  res == -2)  {  item = output;  }

    /** We can't have a non branching node in the wanted direction => iteration is finished. */
    if (res <= 0)  {  isDone = true;  }

}};

template<typename Node, typename Edge, typename GraphDataVariant> 
GraphIterator<Edge> GraphTemplate<Node, Edge, GraphDataVariant>::getSimpleEdgeIterator (Node& node, Direction dir) const
{
    return GraphIterator<Edge> (new EdgeSimplePathIterator<Node, Edge, Functor_getSimpleEdgeIterator<Node, Edge, GraphDataVariant>, GraphDataVariant>(*this, node, dir, Functor_getSimpleEdgeIterator<Node, Edge, GraphDataVariant>()));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant> 
std::set<BranchingNode_t<Node> > GraphTemplate<Node, Edge, GraphDataVariant>::neighbors (typename std::set<BranchingNode_t<Node> >::iterator first, typename std::set<BranchingNode_t<Node> >::iterator last) const
{
#if 0
    std::set<BranchingNode> result;
    for (auto it=first; it!=last; ++it)
    {
        GraphVector<BranchingNode> neighbors = this->neighborsBranching(it->kmer);
        for (size_t i=0; i<neighbors.size(); i++)  { result.insert (neighbors[i]); }
    }
    return result;
#else

    static const size_t nbThread = 8;
    static const size_t nbThreshold = nbThread*1;

    std::set<BranchingNode_t<Node> > result;

    size_t nb = std::distance (first, last);

    if (nb >= nbThreshold)
    {
        typename std::set<BranchingNode_t<Node> >::iterator begin = first;
        typename std::set<BranchingNode_t<Node> >::iterator end;


        int nbPerThread = nb/nbThread;

        vector<pair<typename std::set<BranchingNode_t<Node> >::iterator, typename std::set<BranchingNode_t<Node> >::iterator> > iteratorPairs;

        class Cmd : public tools::dp::ICommand, public system::SmartPointer
        {
        public:
            Cmd (const GraphTemplate<Node, Edge, GraphDataVariant>& graph, const pair<typename std::set<BranchingNode_t<Node> >::iterator, typename std::set<BranchingNode_t<Node> >::iterator>& range)
                : graph(graph), range(range)
            {
                result.reserve (std::distance(range.first,range.second)*8);
            }

            ~Cmd () {}

            void execute ()
            {
                for (typename std::set<BranchingNode_t<Node> >::iterator it=range.first; it!=range.second; ++it)
                {
                    GraphVector<BranchingNode_t<Node> > neighbors = graph.neighborsBranching(it->kmer);
                    for (size_t i=0; i<neighbors.size(); i++)  { result.push_back (neighbors[i]); }
                }
            }

            vector<BranchingNode_t<Node> >& get() { return result; }

        private:
            const GraphTemplate<Node, Edge, GraphDataVariant>& graph;
            pair<typename std::set<BranchingNode_t<Node> >::iterator, typename std::set<BranchingNode_t<Node> >::iterator> range;
            vector<BranchingNode_t<Node> > result;
        };

        vector<tools::dp::ICommand*> cmds;

        while (end != last)
        {
            end = begin;
            std::advance (end, std::min (nbPerThread, (int)std::distance(end,last)));
            iteratorPairs.push_back (std::make_pair(begin, end) );

			tools::dp::ICommand* cmd = new Cmd (*this, std::make_pair(begin, end));
			cmd->use();
			cmds.push_back (cmd);
            begin = end;
        }

        tools::dp::impl::Dispatcher().dispatchCommands(cmds, 0);

        for (size_t i=0; i<cmds.size(); i++)
        {
            vector<BranchingNode_t<Node> >& current = ((Cmd*)cmds[i])->get();

            result.insert (current.begin(), current.end());
            cmds[i]->forget();
        }
    }
    else
    {
        for (typename std::set<BranchingNode_t<Node> >::iterator it=first; it!=last; ++it)
        {
            GraphVector<BranchingNode_t<Node> > neighbors = this->neighborsBranching(it->kmer);
            for (size_t i=0; i<neighbors.size(); i++)  { result.insert (neighbors[i]); }
        }
    }

    return result;
#endif
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : what's this? is it used? no documentation
*********************************************************************/
#if 0 // this code crashes clang, so i disable it for now; actually, i still dont know whether this code is used anywhere
template<typename Node, typename Edge, typename GraphDataVariant> 
struct mutate_visitor : public boost::static_visitor<GraphVector<Node> >    {

    const Node& node;  size_t idx;  int mode;

    mutate_visitor (const Node& node, size_t idx, int mode) : node(node), idx(idx), mode(mode){}

    template<size_t span>  GraphVector<Node> operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::Type           Type;
        typedef typename Kmer<span>::ModelCanonical Model;

        GraphVector<Node> result;
        size_t nbMutations = 0;

        size_t kmerSize = data._model->getKmerSize();

        Model* model = data._model;

        Type kmer = node.template getKmer<Type>();

        size_t idxReverse = kmerSize - 1 - idx;

        if (node.strand == STRAND_FORWARD)
        {
            Nucleotide ntInit = (Nucleotide) (node.kmer[idxReverse]);
            size_t     nt0    = mode==0 ? 0 : ntInit+1;

            Type resetMask;
            resetMask.setVal(3);
            resetMask <<= 2*idxReverse;
            resetMask = ~resetMask;


            for (size_t nt=nt0; nt<4; nt++)
            {
                Type single_nt;
                single_nt.setVal(nt);
                single_nt <<= 2*idxReverse;
 
                Type direct = (kmer & resetMask) | single_nt;
                Type rev    = model->reverse(direct);

                if (direct < rev)
                {
                    if (data.contains(direct))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = direct;
                        node.strand = STRAND_FORWARD;
                        node.mphfIndex = 0;
                        // would need to add iterationRank if it's enabled 
                    }
                }
                else
                {
                    if (data.contains(rev))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = rev;
                        node.strand = STRAND_REVCOMP;
                        node.mphfIndex = 0;
                        // would need to add iterationRank if it's enabled 
                    }
                }
            } /* end of or (size_t nt=nt0; */
        }
        else
        {
            Nucleotide ntInit = reverse ((Nucleotide) (node.kmer[idx]));
            size_t nt0 = mode==0 ? 0 : ntInit+1;

            Type resetMask;
            resetMask.setVal(3);
            resetMask <<= 2*idx;
            resetMask = ~resetMask;

            for (size_t nt=nt0; nt<4; nt++)
            {
                Type single_rev_nt;
                single_rev_nt.setVal(reverse((Nucleotide)nt));
                single_rev_nt <<= 2*idx;
                Type direct = (kmer & resetMask) +  single_rev_nt;
                Type rev    = model->reverse(direct);

                if (direct < rev)
                {
                    if (data.contains(direct))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = direct;
                        node.strand = STRAND_REVCOMP;
                        node.mphfIndex = 0;
                        // would need to add iterationRank if it's enabled 
                    }
                }
                else
                {
                    if (data.contains(rev))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = rev;
                        node.strand = STRAND_FORWARD;
                        node.mphfIndex = 0;
                        // would need to add iterationRank if it's enabled 
                    }
                }
            }
        }

        result.resize(nbMutations);

        return result;
    }
};

/** */
template<typename Node, typename Edge, typename GraphDataVariant> 
GraphVector<Node> GraphTemplate<Node, Edge, GraphDataVariant>::mutate (const Node& node, size_t idx, int mode) const
{
    return boost::apply_visitor (mutate_visitor<Node, Edge, GraphDataVariant>(node,idx,mode),  *(GraphDataVariant*)_variant);
}
#endif

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Node, typename Edge, typename GraphDataVariant> 
struct getNT_visitor : public boost::static_visitor<Nucleotide>    {

    const Node& node;  size_t idx;

    getNT_visitor (const Node& node, size_t idx) : node(node), idx(idx) {}

    template<size_t span>  Nucleotide operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        size_t kmerSize = data._model->getKmerSize();

        if (node.strand == STRAND_FORWARD)  { return (Nucleotide) (node.kmer[kmerSize-1-idx]); }
        else                                { return reverse ((Nucleotide) (node.kmer[idx])); }
    }
};

/** */
template<typename Node, typename Edge, typename GraphDataVariant> 
Nucleotide GraphTemplate<Node, Edge, GraphDataVariant>::getNT (const Node& node, size_t idx) const
{
    return boost::apply_visitor (getNT_visitor<Node, Edge, GraphDataVariant>(node,idx),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/

template<size_t span, typename Node_in>
unsigned long getNodeIndex (const GraphData<span>& data, Node_in& node)
{
    typedef typename Kmer<span>::Type  Type;
    if (node.mphfIndex != 0) // well, this means the 0th node mphf index isn't cached. good enough for me.
        return node.mphfIndex;

    //if (data._abundance == NULL)
    //    throw system::Exception ("Error! getNodeIndex called but MPHF not constructed");
    /* I'm hesitant to put a check (possible branch) in such a critical code. let's check at a higher level */

    /** We get the specific typed value from the generic typed value. */
    Type value = node.template getKmer<Type>();

    // this code was used to make sure that getNodeIndex was always called on a canonical kmer. it is.
#if 0
    {
        size_t      kmerSize = data._model->getKmerSize();
        Type value2 =  revcomp (value, kmerSize);
        if (value2 < value)
        {
            std::cout<<"Error: getNodeIndex has been called on a node that has a kmer that's not canonical"<< std::endl;
            exit(1);
        }
    }
#endif

    // we use _abundance as the mphf. we could also use _nodestate but it might be null if disabled by disableNodeState()
    unsigned long hashIndex = (*(data._abundance)).getCode(value);

    node.mphfIndex = hashIndex;
    
    //std::cout<<"getNodeIndex : " << value <<  " : " << hashIndex << std::endl;
    return hashIndex;
}


/* this whole visitor pattern thing in the GraphTemplate<Node, Edge, GraphDataVariant>..
  it is used to support querying the right graph variant (the one that corresponds to the adequate kmer size)
*/
template<typename Node, typename Edge, typename GraphDataVariant> 
struct queryAbundance_visitor : public boost::static_visitor<int>    {

    Node& node;

    queryAbundance_visitor (Node& node) : node(node){}

    template<size_t span>  int operator() (const GraphData<span>& data) const
    {
        unsigned long hashIndex = getNodeIndex<span>(data, node);
    	if(hashIndex == ULLONG_MAX) return 0; // node was not found in the mphf 

        int value = data._abundance->abundanceAt(hashIndex); // uses discretized abundance

        return value;
    }
};

/** */
template<typename Node, typename Edge, typename GraphDataVariant> 
int GraphTemplate<Node, Edge, GraphDataVariant>::queryAbundance (Node& node) const
{
    return boost::apply_visitor (queryAbundance_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}

/* 
/ a node state, using the MPHF, is either:
/ 0: unmarked (normal state)
/ 1: marked (already in an output unitig/contig)
/ 2: deleted (in minia: tips and collapsed bubble paths will be deleted)
/ 3: complex (branching or deadend, unmarked)
*/

template<typename Node, typename Edge, typename GraphDataVariant> 
struct queryNodeState_visitor : public boost::static_visitor<int>    {

    Node& node;

    queryNodeState_visitor (Node& node) : node(node){}

    template<size_t span>  int operator() (const GraphData<span>& data)  const
    {
        unsigned long hashIndex = getNodeIndex<span>(data, node);
    	if(hashIndex == ULLONG_MAX) return 0; // node was not found in the mphf 

        unsigned char value = (*(data._nodestate)).at(hashIndex / 2);

        if (hashIndex % 2 == 1)
            value >>= 4;

        value &= 0xF;

        return value;
    }
};

/** */
template<typename Node, typename Edge, typename GraphDataVariant>
int GraphTemplate<Node, Edge, GraphDataVariant>::queryNodeState (Node& node) const 
{
    return boost::apply_visitor (queryNodeState_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}



template<typename Node, typename Edge, typename GraphDataVariant> 
struct setNodeState_visitor : public boost::static_visitor<int>    {

    Node& node;

    int state;

    setNodeState_visitor (Node& node, int state) : node(node), state(state) {}

    template<size_t span> int operator() (const GraphData<span>& data)  const
    {
        unsigned long hashIndex = getNodeIndex<span>(data, node);
    	if(hashIndex == ULLONG_MAX) return 0; // node was not found in the mphf 

        unsigned char &value = (*(data._nodestate)).at(hashIndex / 2);

        int maskedState = state & 0xF;

        if (hashIndex % 2 == 1)
        {
            value &= 0xF;
            value |= (maskedState << 4);
        }
        else
        {
            value &= 0xF0;
            value |= maskedState;
        }

        return 0;
    }
};

/** */
template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::setNodeState (Node& node, int state) const 
{
    boost::apply_visitor (setNodeState_visitor<Node, Edge, GraphDataVariant>(node, state),  *(GraphDataVariant*)_variant);
}

template<typename Node, typename Edge, typename GraphDataVariant>
struct resetNodeState_visitor : public boost::static_visitor<int>    {

    resetNodeState_visitor () {}

    template<size_t span> int operator() (const GraphData<span>& data) const
    {
        (*(data._nodestate)).clearData();
        return 0;
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::resetNodeState() const
{
    boost::apply_visitor (resetNodeState_visitor<Node, Edge, GraphDataVariant>(),  *(GraphDataVariant*)_variant);
}


template<typename Node, typename Edge, typename GraphDataVariant>
struct disableNodeState_visitor : public boost::static_visitor<int>    {

    disableNodeState_visitor () {}

    template<size_t span> int operator() (GraphData<span>& data) const
    {
        data._nodestate = NULL;
        return 0;
    }
};

/* this can be useful for benchmarking purpose.
 * because contains() always does a MPHF query to check if the node is deleted or not
 * unless _nodestate is NULL.
 * thus, this functions sets _nodestate to NULL.
 * NOTE: irreversible!
 */
template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::disableNodeState() const
{
    boost::apply_visitor (disableNodeState_visitor<Node, Edge, GraphDataVariant>(),  *(GraphDataVariant*)_variant);
}

template<typename Node, typename Edge, typename GraphDataVariant> 
bool GraphTemplate<Node, Edge, GraphDataVariant>::isNodeDeleted(Node& node) const
{
    return (!checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE)) || (((queryNodeState(node) >> 1) & 1) == 1);
}


// direct access to MPHF index

template<typename Node, typename Edge, typename GraphDataVariant> 
struct nodeMPHFIndex_visitor : public boost::static_visitor<unsigned long>    {

    Node& node;

    nodeMPHFIndex_visitor (Node& node) : node(node) {}

    template<size_t span> unsigned long operator() (const GraphData<span>& data)  const
    {
        return getNodeIndex<span>(data, node);
    }
};

template<typename Node, typename Edge, typename GraphDataVariant> 
unsigned long GraphTemplate<Node, Edge, GraphDataVariant>::nodeMPHFIndex(Node& node) const 
{
    if (!checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE))
       return 0;
    return boost::apply_visitor (nodeMPHFIndex_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}

/* debug function, only for profiling */
template<typename Node, typename Edge, typename GraphDataVariant> 
struct nodeMPHFIndex_visitorDummy : public boost::static_visitor<unsigned long>    {

    const Node& node;

    nodeMPHFIndex_visitorDummy (const Node& node) : node(node) {}

    template<size_t span> unsigned long operator() (const GraphData<span>& data) const
    {
        /** We get the specific typed value from the generic typed value. */
        unsigned long hashIndex = 0;
        return hashIndex;
    }
};

template<typename Node, typename Edge, typename GraphDataVariant> 
unsigned long GraphTemplate<Node, Edge, GraphDataVariant>::nodeMPHFIndexDummy(Node& node) const /* debug function, only for profiling*/
{
    if (!checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE))
       return 0;
    return boost::apply_visitor (nodeMPHFIndex_visitorDummy<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}


/*****************************************************************************
 *
 * graph adjacency functions 
 *
 *
 * */


template<typename Node, typename Edge, typename GraphDataVariant> 
struct getAdjacency_visitor : public boost::static_visitor<unsigned char&>    {

    Node& node;

    getAdjacency_visitor (Node& node) : node(node) {}

    template<size_t span> unsigned char& operator() (const GraphData<span>& data) const
    {
        unsigned long hashIndex= getNodeIndex<span>(data, node);
    	if(hashIndex == ULLONG_MAX) 
        { // node was not found in the mphf: complain a return a dummy value
            std::cout << "getAdjacency called for node not in MPHF" << std::endl; 
            return (*(data._adjacency)).at(0);
        }

        unsigned char &value = (*(data._adjacency)).at(hashIndex);
        //std::cout << "hashIndex " << hashIndex << " value " << (int)value << std::endl;;

        return value;
    }
};

/* initiate the adjacency map */
template<typename Node, typename Edge, typename GraphDataVariant> 
struct allocateAdjacency_visitor : public boost::static_visitor<void>    {

    template<size_t span> void operator() (const GraphData<span>& data) const
    {
        data._adjacency->useHashFrom(data._abundance); // use abundancemap's MPHF, and allocate 8 bits per element for the adjacency map
    }
};


/* precompute the graph adjacency information using the MPHF
 * this should be much faster than querying the bloom filter
 * also, maybe one day, it will replace it
 */
template<typename Node, typename Edge, typename GraphDataVariant> 
void GraphTemplate<Node, Edge, GraphDataVariant>::precomputeAdjacency(unsigned int nbCores, bool verbose) 
{
    ProgressGraphIteratorTemplate<Node, ProgressTimerAndSystem> itNode (iterator(), "precomputing adjacency", verbose);
    
    bool hasMPHF = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE;
    if (!hasMPHF)
    {
        throw system::Exception ("Cannot precompute adjacency information - MPHF was not constructed");
        return;
    }

    /* allocate the adjacency map */
    boost::apply_visitor (allocateAdjacency_visitor<Node, Edge, GraphDataVariant>(),  *(GraphDataVariant*)_variant);

    Dispatcher dispatcher (nbCores); 

    // nt2int is defined as static in Graph.hpp
    nt2bit[NUCL_A] = 1;
    nt2bit[NUCL_C] = 2;
    nt2bit[NUCL_T] = 4;
    nt2bit[NUCL_G] = 8;

    dispatcher.iterate (itNode, [&] (Node& node)        {

            unsigned char &value = boost::apply_visitor (getAdjacency_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
            value = 0;

            // in both directions
            for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
            {
                GraphVector<Edge> neighbors = this->neighborsEdge(node, dir);
                for (unsigned int i = 0; i < neighbors.size(); i++)
                {
                    u_int8_t bit = nt2bit[neighbors[i].nt];
                    value |= bit << (dir == DIR_INCOMING ? 4 : 0);
                    //std::cout << "setting bit " << (int)bit << " shifted " << ((int)bit << (dir == DIR_INCOMING ? 4 : 0)) << " for nt " << (int)(neighbors[i].nt) << std::endl;
                }


                //std::cout << "node " << this->toString(node) << " has " << neighbors.size() << " neighbors in direction " << (dir == DIR_INCOMING ? "incoming" : "outcoming") << " value is now " << (int)value <<  std::endl;
            }
            
    }); // end of parallel node iterate

    setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE);
    
    // do a sanity check, to see if adjacency information matches bloom information 
    bool adjSanityCheck = false;
    if (adjSanityCheck)
    {
        ProgressGraphIteratorTemplate<Node, ProgressTimerAndSystem> itNode2 (iterator(), "checking adjacency", verbose);
        dispatcher.iterate (itNode2, [&] (Node& node)        {
                // in both directions
                for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
                {
                    for (int otherStrand = 0; otherStrand <= 2; otherStrand++) // test both strands
                    {
                        node.strand = (otherStrand == 0) ? STRAND_FORWARD: STRAND_REVCOMP;
                        if (this->debugCompareNeighborhoods(node, dir, "bad adjacency check")) exit(1);
                    }
                }
        }); // end of parallel node iterate
    }


    // TODO delete _container here
}

// now deleteNode depends on getNodeAdjacency
template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::deleteNode (Node& node) const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    if (hasAdjacency)
    {
        /* absolutely need to update adjacency information of neighboring nodes */
        // todo: this can be made more efficient, but I wanted quick bug-free code
        GraphVector<Edge> neighbs = this->neighborsEdge(node);

        // sanity checks regarding adjacency information
        //if (debugCompareNeighborhoods(node, DIR_INCOMING, "node deletion, incoming neighbors")) exit(1);
        //if (debugCompareNeighborhoods(node, DIR_OUTCOMING, "node deletion, outcoming neighbors")) exit(1);

        for (size_t i = 0; i < neighbs.size(); i++)
        {
            Node neighbor = neighbs[i].to;

            bool deleted = false;

            // also, since in precomputeAdjacency i'm using DIR_INCOMING and DIR_OUTCOMING, let's do the same here. 
            // (gatb seems to ahve different edge creation behavior when neighbors(node) and neighbors(node,DIR_INCOMING) is used)
            for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
            {
                GraphVector<Edge> neighbs_of_neighbs = this->neighborsEdge(neighbor, dir);
                
                // sanity check 
                //if (debugCompareNeighborhoods(neighbor,  dir, "node deletion, neighbs_of_neighbs"))
                //{
                //    std::cout << "some debug info, node " << this->toString(node) << ", neighbor " << this->toString(neighbor) << " strand: " << ((neighbor.strand == STRAND_FORWARD)? "forward":"reverse") << " direction: " <<  (dir == DIR_INCOMING ? "incoming": "outcoming") << std::endl; exit(1);
                //}

                for (size_t j = 0; j  < neighbs_of_neighbs.size(); j++)
                {   
                    Nucleotide nt = neighbs_of_neighbs[j].nt;
                    Node neigh_of_neigh = neighbs_of_neighbs[j].to;
                    
                    if (neigh_of_neigh == node)
                    {
                        unsigned char& value = boost::apply_visitor (getAdjacency_visitor<Node, Edge, GraphDataVariant>(neighbor),  *(GraphDataVariant*)_variant);
                        u_int8_t bit = nt2bit[nt];
                            
                        bool forwardStrand = (neighbor.strand == STRAND_FORWARD);
                        int shift;
                        if (forwardStrand) 
                        {
                            shift = (dir == DIR_INCOMING) ? 4 : 0;
                        }
                        else
                        {
                            shift = (dir == DIR_OUTCOMING) ? 4 : 0;
                            bit = ((bit & 3) << 2) | ((bit >> 2) & 3);
                        }
                        //std::cout << "deleting node with value " << (int)value << " and dir :" << (dir == DIR_INCOMING ? "incoming": "outcoming") << ": neighbor" << ((neighbor.strand==STRAND_REVCOMP) ? "(r)":"")<<" " << this->toString(neighbor) << " --(nt=" << nt << ")--> neigh_of_neigh"  << ((neigh_of_neigh.strand==STRAND_REVCOMP) ? "(r)":"")<< " " << this->toString(neigh_of_neigh) << std::endl;

                        if (((value >> shift) & bit) == 0) // TODO remove this check if no problem after a while
                        {
                            std::cout << "Error while deleting node " <<  this->toString(node) << ": neighbor" << ((neighbor.strand==STRAND_REVCOMP) ? "(r)":"")<<" " << this->toString(neighbor) << " --(nt=" << nt << ")--> neigh_of_neigh"  << ((neigh_of_neigh.strand==STRAND_REVCOMP) ? "(r)":"")<< " " << this->toString(neigh_of_neigh) << " and dir :" << (dir == DIR_INCOMING ? "incoming": "outcoming") << ", value " << (int)value << std::endl;
                            exit(1);
                        }
                        value ^= (bit << shift);
                        
                        deleted = true;
                    }
                }
            }
            if (deleted == false)
            {
                std::cout << "[bug detected in adjacency data structure] was supposed to delete node "<< this->toString(node) <<" but couldn't find it :(." << std::endl;
                std::cout << "assembly may be inaccurate around this kmer. will continue assembling anyway" << std::endl;
            }

        }
    }

    // a little sanitycheck
#if 0
    Vector<Edge> neighbs = this->neighborsEdge(node);
#endif

    if (checkState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_MPHF_DONE))
        setNodeState(node, 2);

    // another test
#if 0
    for (size_t i = 0; i < neighbs.size(); i++)
    {
        Node neighbor = neighbs[i].to;
        for (Direction dir=DIR_OUTCOMING; dir<DIR_END; dir = (Direction)((int)dir + 1) )
            if (debugCompareNeighborhoods(neighbor,dir,"post delete node")) exit(1);
    }
#endif

    // update cached branching nodes information now
    bool _cacheNonSimpleNodes = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_NONSIMPLE_CACHE;
    if (_cacheNonSimpleNodes)
    {
        cacheNonSimpleNodeDelete(node); // so in case of a tip, will delete the tip and add the next kmer as non-branching, which will be in turn deleted. not that efficient, but will do for now.
        GraphVector<Edge> neighbs = this->neighborsEdge(node);
        for (size_t i = 0; i < neighbs.size(); i++)
        {
            Node neighbor = neighbs[i].to;
            if (isBranching(neighbor))
                cacheNonSimpleNode(neighbor);
        }
    }
}

template<typename Node, typename Edge, typename GraphDataVariant>
bool GraphTemplate<Node, Edge, GraphDataVariant>::debugCompareNeighborhoods(Node& node, Direction dir, std::string prefix) const
{
    bool hasAdjacency = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_ADJACENCY_DONE;
    if (!hasAdjacency) return false;

    GraphVector<Edge> neighborsAdj = boost::apply_visitor (getItems_visitor<Node, Edge, Edge, Functor_getEdges<Node, Edge, GraphDataVariant>, GraphDataVariant>(node, dir, true , Functor_getEdges<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant);
    GraphVector<Edge> neighborsBloom = boost::apply_visitor (getItems_visitor<Node, Edge, Edge, Functor_getEdges<Node, Edge, GraphDataVariant>, GraphDataVariant>(node, dir, false, Functor_getEdges<Node, Edge, GraphDataVariant>()),  *(GraphDataVariant*)_variant); // without adjacency (copied neighborsEdge code because hasAdjacency parameter isn't exposed)

    if (neighborsAdj.size() != neighborsBloom.size())
    {
        std::cout << prefix << " #### mismatched number of neighbors: node " << this->toString(node) <<  " dir " << (dir == DIR_OUTCOMING? "outcoming" :" incoming") << ", " << neighborsAdj.size() << " neighborsAdj, " << neighborsBloom.size()  << " neighborsBloom" << std::endl;
        return true;
    }

    bool bad=false;
    for (unsigned int j = 0; j < neighborsAdj.size(); j++)
    {
        if (neighborsAdj[j].to != neighborsBloom[j].to)
        {
            std::cout << prefix << "#### node " << this->toString(node) << " neighbor error (out of " << neighborsAdj.size() << " neighbor(s)): neighborsAdj[j].to = " << this->toString(neighborsAdj[j].to) << ", neighborsBloom[j].to = " << this->toString(neighborsBloom[j].to) << std::endl;
            bad = true;
        }

        if (neighborsAdj[j].from != neighborsBloom[j].from)
        {
            std::cout << prefix << "#### node " << this->toString(node) << " neighbor error: neighborsAdj[j].from = " << this->toString(neighborsAdj[j].from) << ", neighborsBloom[j].from = " << this->toString(neighborsBloom[j].from) << std::endl;
            bad = true;
        }
        if (neighborsAdj[j].nt != neighborsBloom[j].nt)
        {
            std::cout << prefix << "#### node " << this->toString(node) << " nucleotide error: neighborsAdj[j].nt=" << neighborsAdj[j].nt << ", neighborsBloom[j].nt=" << neighborsBloom[j].nt << std::endl;
            bad = true;
        }

    }
    return bad;
}

// TODO: it makes sense someday to introduce a graph._nbCore parameter, because this function, simplify() and precomputeAdjacency() all want it
template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::deleteNodesByIndex(vector<bool> &bitmap, int nbCores, gatb::core::system::ISynchronizer* synchro) const
{
    GraphIterator<Node> itNode = this->iterator();
    Dispatcher dispatcher (nbCores); 

    dispatcher.iterate (itNode, [&] (Node& node)        {

        unsigned long i = this->nodeMPHFIndex(node); 

        if (bitmap[i])
        {
            // make this deletion atomic with respect to other node deletions
            if (synchro)
                synchro->lock();

            this->deleteNode(node);

            if (synchro)
                synchro->unlock();
        }
    }); // end of parallel node iteration
}

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::simplify(unsigned int nbCores, bool verbose)
{
        Simplifications<GraphTemplate<Node, Edge, GraphDataVariant>,Node,Edge> 
            graphSimplifications(this, nbCores, verbose);
        graphSimplifications.simplify();
}



template<typename Node, typename Edge, typename GraphDataVariant> 
struct cacheNonSimpleNode_visitor : public boost::static_visitor<void>    {
    const Node& node;

    cacheNonSimpleNode_visitor (const Node& node) : node(node) {}

    template<size_t span> void operator() (const GraphData<span>& data) const
    {
        typedef typename Kmer<span>::Type Type;
        Type value = node.template getKmer<Type>();
        if (data._nodecache->find(value) != data._nodecache->end())
            return; // don't overwrite existing cached node
        // add dummy cached sequence info
        std::pair<char,string> dummy (-1,"");
        (*(data._nodecache))[value] = dummy;
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::cacheNonSimpleNode(const Node& node) const 
{
    bool _cacheNonSimpleNodes = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_NONSIMPLE_CACHE;
    if (!_cacheNonSimpleNodes)
        return; // don't do anything if we don't cache nodes
    boost::apply_visitor (cacheNonSimpleNode_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}

template<typename Node, typename Edge, typename GraphDataVariant> 
struct cacheNonSimpleNodeDelete_visitor : public boost::static_visitor<void>    {

    const Node& node;

    cacheNonSimpleNodeDelete_visitor (const Node& node) : node(node) {}

    template<size_t span> void operator() (const GraphData<span>& data) const
    {
        typedef typename Kmer<span>::Type Type;
        Type value = node.template getKmer<Type>();
        typename GraphData<span>::NodeCacheMap::iterator pos = data._nodecache->find(value);
        if (pos != data._nodecache->end())
            data._nodecache->erase(pos);
      //  else
      //      std::cout << "Warning: attempting to delete a node from the cache, but it wasn't present." << std::endl;
      /* not printing that warning because will call the deleter on any node that's deleted, not just branching ones */
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::cacheNonSimpleNodeDelete(const Node& node) const 
{
    bool _cacheNonSimpleNodes = getState() & GraphTemplate<Node, Edge, GraphDataVariant>::STATE_NONSIMPLE_CACHE;
    if (!_cacheNonSimpleNodes)
        return; // don't do anything if we don't cache nodes
    boost::apply_visitor (cacheNonSimpleNodeDelete_visitor<Node, Edge, GraphDataVariant>(node),  *(GraphDataVariant*)_variant);
}

template<typename Node, typename Edge, typename GraphDataVariant> 
struct allocateNonSimpleNodeCache_visitor : public boost::static_visitor<void>    {

    template<size_t span> void operator() (GraphData<span>& data) const 
    {
        data.setNodeCache(new typename GraphData<span>::NodeCacheMap);
    }
};


template<typename Node, typename Edge, typename GraphDataVariant>
void GraphTemplate<Node, Edge, GraphDataVariant>::cacheNonSimpleNodes(unsigned int nbCores, bool verbose) 
{
    boost::apply_visitor (allocateNonSimpleNodeCache_visitor<Node, Edge, GraphDataVariant>(),  *(GraphDataVariant*)_variant);
    setState(GraphTemplate<Node, Edge, GraphDataVariant>::STATE_NONSIMPLE_CACHE);
    GraphIterator<Node> itNode = this->iterator();
    Dispatcher dispatcher (nbCores); 
    system::ISynchronizer* synchro = system::impl::System::thread().newSynchronizer();
    unsigned long nbCachedNodes = 0;
    dispatcher.iterate (itNode, [&] (Node& node)        {
        if (this->isNodeDeleted(node)) return; // test
        if (this->isBranching(node))
        {
            synchro->lock();
            this->cacheNonSimpleNode(node);
            __sync_fetch_and_add(&nbCachedNodes,1);
            synchro->unlock();
        }
    }); // end of parallel node iteration
    std::cout << "Cached " << nbCachedNodes << " non-simple nodes" << std::endl;
}

template<typename Node, typename Edge, typename GraphDataVariant>
struct cached_nodes_visitor : public boost::static_visitor<tools::dp::ISmartIterator<Node>*>
{
    const GraphTemplate<Node, Edge, GraphDataVariant>& graph;
    cached_nodes_visitor (const GraphTemplate<Node, Edge, GraphDataVariant>& graph) : graph(graph) {}

    // we should really use STL iterators in the next rewrite.
    template<size_t span>  tools::dp::ISmartIterator<Node>* operator() (const GraphData<span>& data) const
    {
        class CachedNodeIterator : public tools::dp::ISmartIterator<Node>
        {
        public:
            CachedNodeIterator (typename GraphData<span>::NodeCacheMap *nodecache)
                : _nodecache(nodecache), _rank(0), _isDone(true) {  
                    this->_item->strand = STRAND_FORWARD;  // iterated nodes are always in forward strand.
                }

            u_int64_t rank () const { return _rank; }

            /** \copydoc  Iterator::first */
            void first()
            {
                _it = _nodecache->begin();
                _rank   = 0;
                _isDone = _it == _nodecache->end();

                if (!_isDone)
                {
                    this->_rank ++;
                    this->_item->kmer      = _it->first;
                    this->_item->abundance = 0; // not recorded
                    this->_item->mphfIndex = 0;
                    this->_item->iterationRank = this->_rank;
                }
            }

            /** \copydoc  Iterator::next */
            void next()
            {
                _it++;
                _isDone = _it == _nodecache->end();
                if (!_isDone)
                {
                    // NOTE: doesn't check if node is deleted (as it would be expensive to compute MPHF index)
                    this->_rank ++;
                    this->_item->kmer      = _it->first;
                    this->_item->abundance = 0; // not recorded
                    this->_item->mphfIndex = 0;
                    this->_item->iterationRank = this->_rank;
                }
            }

            /** \copydoc  Iterator::isDone */
            bool isDone() { return _isDone;  }

            /** \copydoc  Iterator::item */
            Node& item ()  {  return *(this->_item);  }

            /** */
            u_int64_t size () const { return _nodecache->size(); }

        private:
            typename GraphData<span>::NodeCacheMap *_nodecache;
            typename GraphData<span>::NodeCacheMap::iterator _it;

            u_int64_t _rank;
            bool      _isDone;
            u_int64_t _nbItems;
        };

        return new CachedNodeIterator (data._nodecache);
    }
};

template<typename Node, typename Edge, typename GraphDataVariant>
GraphIterator<Node> GraphTemplate<Node, Edge, GraphDataVariant>::iteratorCachedNodes() const
{
    return GraphIterator<Node> (boost::apply_visitor (cached_nodes_visitor<Node,Edge, GraphDataVariant>(*this),  *(GraphDataVariant*)_variant));
}


// instantiation
// uses Node and Edge as defined in Graph.hpp (legacy GATB compatibility, when Graph was not templated)
template class GraphTemplate<Node, Edge, GraphDataVariant>; 

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
