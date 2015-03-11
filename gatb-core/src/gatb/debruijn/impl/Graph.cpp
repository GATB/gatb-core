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

#include <gatb/tools/collections/impl/ContainerSet.hpp>

#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/LibraryInfo.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/tools/misc/impl/Tool.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/BloomAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithmFactory.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/kmer/impl/MPHFAlgorithm.hpp>
#include <gatb/tools/collections/impl/MPHF.hpp>

using namespace std;

using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

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

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

/* We define a structure that holds all the necessary stuff for implementing the graph API.
 *  Here, the structure is templated by the span (ie. the kmer max size).
 *
 *  This structure is the basis for defining a boost::variant with all required span
 *  template instantiations.
 */
template<size_t span>
struct GraphData
{
    /** Shortcuts. */
    typedef typename Kmer<span>::ModelCanonical Model;
    typedef typename Kmer<span>::Type           Type;
    typedef typename Kmer<span>::Count          Count;
    typedef typename MPHFAlgorithm<span>::Map   Map;

    /** Constructor. */
    GraphData () : _model(0), _solid(0), _container(0), _branching(0), _abundance(0) {}

    /** Destructor. */
    ~GraphData ()
    {
        setModel     (0);
        setSolid     (0);
        setContainer (0);
        setBranching (0);
        setAbundance (0);
    }

    /** Constructor (copy). */
    GraphData (const GraphData& d) : _model(0), _solid(0), _container(0), _branching(0), _abundance(0)
    {
        setModel     (d._model);
        setSolid     (d._solid);
        setContainer (d._container);
        setBranching (d._branching);
        setAbundance (d._abundance);
    }

    /** Assignment operator. */
    GraphData& operator= (const GraphData& d)
    {
        if (this != &d)
        {
            setModel     (d._model);
            setSolid     (d._solid);
            setContainer (d._container);
            setBranching (d._branching);
            setAbundance (d._abundance);
        }
        return *this;
    }

    /** Required attributes. */
    Model*                _model;
    Partition<Count>*     _solid;
    IContainerNode<Type>* _container;
    Collection<Count>*    _branching;
    Map*                  _abundance;

    /** Setters. */
    void setModel       (Model*                 model)      { SP_SETATTR (model);     }
    void setSolid       (Partition<Count>*      solid)      { SP_SETATTR (solid);     }
    void setContainer   (IContainerNode<Type>*  container)  { SP_SETATTR (container); }
    void setBranching   (Collection<Count>*     branching)  { SP_SETATTR (branching); }
    void setAbundance   (Map*                   abundance)  { SP_SETATTR (abundance); }

    /** Shortcut. */
    bool contains (const Type& item)  const  {  return _container->contains (item);  }
};

/********************************************************************************/

/* This definition is the basis for having a "generic" Graph class, ie. not relying on a template
 * parameter.
 *
 * This is done through a boost::variant; actually, we use a limited number of variant, corresponding
 * to some maximum kmer sizes.
 */
typedef boost::variant <
    GraphData<KSIZE_1>,
    GraphData<KSIZE_2>,
    GraphData<KSIZE_3>,
    GraphData<KSIZE_4>
>  GraphDataVariant;

/********************************************************************************/

/** This function set a specific value to the given GraphDataVariant object,
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
static void setVariant (GraphDataVariant& data, size_t kmerSize, size_t integerPrecision=0)
{
	/** We may force the kmer size and not use the optimized KSIZE value. */
	if (integerPrecision > 0)  { kmerSize = integerPrecision*32 - 1; }

    /** Here is the link between the kmer size (or precision) and the specific type to be used for the variant. */
         if (kmerSize < KSIZE_1)  {  data = GraphData<KSIZE_1> (); }
    else if (kmerSize < KSIZE_2)  {  data = GraphData<KSIZE_2> (); }
    else if (kmerSize < KSIZE_3)  {  data = GraphData<KSIZE_3> (); }
    else if (kmerSize < KSIZE_4)  {  data = GraphData<KSIZE_4> (); }
    else { throw system::Exception ("Graph failure because of unhandled kmer size %d", kmerSize); }

    /** We convert the kmer size chosen by the user to the precision value. */
    size_t prec = 1 + kmerSize / 32;
    Integer::setType (prec);
}

/********************************************************************************/

/* This visitor is used to configure a GraphDataVariant object (ie configure its attributes).
 * The information source used to configure the variant is a kmer size and a storage.
 *
 * If the storage is null, only the kmer model is set.
 *
 * If the storage is not null, it is likely coming from a previous graph building (dsk, debloom, ...); we can
 * therefore configure the variant with the items coming from this storage.
 */
struct configure_visitor : public boost::static_visitor<>    {

    const Graph& graph;
    Storage&     storage;

    configure_visitor (const Graph& graph, Storage& storage)  : graph(graph), storage(storage) {}

    template<size_t span>  void operator() (GraphData<span>& data) const
    {
        size_t   kmerSize = graph.getKmerSize();

        /** We create the kmer model. */
        data.setModel (new typename Kmer<span>::ModelCanonical (kmerSize));

        if (graph.getState() & Graph::STATE_SORTING_COUNT_DONE)
        {
            /** We set the iterable for the solid kmers. */
            SortingCountAlgorithm<span> algo (storage);
            graph.getInfo().add (1, algo.getInfo());
            data.setSolid (algo.getSolidCounts());
       }

        if (graph.getState() & Graph::STATE_BLOOM_DONE)
        {
            /** We set the container. */
            BloomAlgorithm<span> algo (storage);
            graph.getInfo().add (1, algo.getInfo());
        }

        if (graph.getState() & Graph::STATE_DEBLOOM_DONE)
        {
            /** We set the container. */
            DebloomAlgorithm<span> algo (storage);
            graph.getInfo().add (1, algo.getInfo());
            data.setContainer (algo.getContainerNode());
        }

        if (graph.getState() & Graph::STATE_BRANCHING_DONE)
        {
            /** We set the branching container. */
            BranchingAlgorithm<span> algo (storage);
            graph.getInfo().add (1, algo.getInfo());
            data.setBranching (algo.getBranchingCollection());
        }

        if ((graph.getState() & Graph::STATE_MPHF_DONE) &&  (graph.getState() & Graph::STATE_SORTING_COUNT_DONE))
        {
            // actually need to get solid kmers (as a way to reconstruct abundance data in the mphf) so i'm calling this as a temporary hack. later: remove that line, just save/load the MPHF abundance array to disk
            SortingCountAlgorithm<span> sortingCount (storage);

            MPHFAlgorithm<span> mphf_algo (
                sortingCount.getStorageGroup(),
                "mphf",
                sortingCount.getSolidCounts(),
                sortingCount.getSolidKmers(),
                false  // build=true, load=false
            );

            data.setAbundance (mphf_algo.getMap());
        }
    }
};

/********************************************************************************/

/* This visitor is used to build a graph. In particular, the data variant of the graph will
 * be configured through this boost visitor.
 *
 * The skeleton of the graph building is the following:
 *  - conversion of the input reads into a binary format
 *  - kmers counting
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
struct build_visitor : public boost::static_visitor<>    {

    Graph& graph; bank::IBank* bank; tools::misc::IProperties* props;

    build_visitor (Graph& aGraph, bank::IBank* aBank, tools::misc::IProperties* aProps)  : graph(aGraph), bank(aBank), props(aProps) {}

    template<size_t span>  void operator() (GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::Type  Type;
        typedef typename Kmer<span>::Count Count;

        LOCAL (bank);

        size_t kmerSize      = props->get(STR_KMER_SIZE)          ? props->getInt(STR_KMER_SIZE)           : 31;
        size_t minimizerSize = props->get(STR_MINIMIZER_SIZE)     ? props->getInt(STR_MINIMIZER_SIZE)      : 8;
        size_t nksMin        = props->get(STR_KMER_ABUNDANCE_MIN) ? props->getInt(STR_KMER_ABUNDANCE_MIN)  : 3;
        size_t nksMax        = props->get(STR_KMER_ABUNDANCE_MAX) ? props->getInt(STR_KMER_ABUNDANCE_MAX)  : 0; // if max<min, we use max=MAX
        size_t minimizerType = props->get(STR_MINIMIZER_TYPE)     ? props->getInt(STR_MINIMIZER_TYPE)      : 0;
        size_t repartitionType = props->get(STR_REPARTITION_TYPE)  ? props->getInt(STR_REPARTITION_TYPE)   : 0;

        string output = props->get(STR_URI_OUTPUT) ?
            props->getStr(STR_URI_OUTPUT)   :
            (props->getStr(STR_URI_OUTPUT_DIR) + "/" + system::impl::System::file().getBaseName (bank->getId()));

        string binaryBankUri =
              System::file().getCurrentDirectory()
            + string("/")
            + System::file().getBaseName (bank->getId())
            + string(".bin");

        DEBUG ((cout << "builGraph for bank '" << bank->getId() << "'"
            << " kmerSize=" << kmerSize
            << " nks=" << nks
            << " output='" << output << "'"
            << endl
        ));

        /** We create the kmer model. */
        data.setModel (new typename Kmer<span>::ModelCanonical (kmerSize));

        /** We add library information. */
        graph.getInfo().add (1, & LibraryInfo::getInfo());

        /************************************************************/
        /*                       Storage creation                   */
        /************************************************************/

        /** We create the storage object for the graph. */
        graph.setStorage (StorageFactory(graph._storageMode).create (output, true, false));

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

        /************************************************************/
        /*                         Sorting count                    */
        /************************************************************/
        KmerSolidityKind solidityKind;  parse (props->getStr(STR_SOLIDITY_KIND), solidityKind);

        /** We create a DSK instance and execute it. */
        SortingCountAlgorithm<span> sortingCount (
            solidStorage,
            bank,
            kmerSize,
            make_pair (nksMin, nksMax),
            props->get(STR_MAX_MEMORY)    ? props->getInt(STR_MAX_MEMORY) : 0,
            props->get(STR_MAX_DISK)      ? props->getInt(STR_MAX_DISK)   : 0,
            props->get(STR_NB_CORES)      ? props->getInt(STR_NB_CORES)   : 0,
            solidityKind,
            props->get(STR_HISTOGRAM_MAX) ? props->getInt(STR_HISTOGRAM_MAX) : 0,
            0,
            minimizerType,
            repartitionType,
            minimizerSize
        );
        executeAlgorithm (sortingCount, *solidStorage, props, graph._info);
        graph.setState(Graph::STATE_SORTING_COUNT_DONE);

        /** We configure the variant. */
        data.setSolid (sortingCount.getSolidCounts());

        /* always print number of solid kmers: this is important information in case a use reports that Graph construction failed/took too long */
        cout << "Found " << sortingCount.getSolidCounts()->getNbItems() << " solid kmers." << endl;

        /** We check that we got solid kmers. */
        if (sortingCount.getSolidCounts()->getNbItems() == 0)  {  return;  /*throw "NO SOLID KMERS FOUND...";*/  }

        /************************************************************/
        /*                         MPHF                             */
        // note: theoretically could be done in parallel to debloom, but both tasks may or may not be IO intensive
        /************************************************************/

        /** We create an instance of the MPHF Algorithm class (why is that a class, and not a function?) and execute it. */
        if (graph._mphfKind != MPHF_NONE)
        {
            MPHFAlgorithm<span> mphf_algo (
                sortingCount.getStorageGroup(),
                "mphf",
                sortingCount.getSolidCounts(),
                sortingCount.getSolidKmers(),
                true  // build=true, load=false
            );
            executeAlgorithm (mphf_algo, graph.getStorage(), props, graph._info);
            data.setAbundance(mphf_algo.getMap());
            graph.setState(Graph::STATE_MPHF_DONE);
        }

        /************************************************************/
        /*                         Bloom                            */
        /************************************************************/
        if (graph.checkState(Graph::STATE_SORTING_COUNT_DONE))
        {
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
                executeAlgorithm (bloomAlgo, graph.getStorage(), props, graph._info);
                graph.setState(Graph::STATE_BLOOM_DONE);
            }
        }

        /************************************************************/
        /*                         Debloom                          */
        /************************************************************/
        if (graph.checkState(Graph::STATE_BLOOM_DONE))
        {
            /** We create a debloom instance and execute it. */
            DebloomAlgorithm<span>* debloom = DebloomAlgorithmFactory<span>::create (
                graph._debloomImpl,
                graph.getStorage(),
                *solidStorage,
                data._solid,
                kmerSize,
                props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 0,
                props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0,
                graph._bloomKind,
                graph._debloomKind
            );
            LOCAL (debloom);

            executeAlgorithm (*debloom, graph.getStorage(), props, graph._info);

            graph.setState(Graph::STATE_DEBLOOM_DONE);

            /** We configure the variant. */
            data.setContainer (debloom->getContainerNode());
        }

        /************************************************************/
        /*                         Branching                        */
        /************************************************************/
        if (graph.checkState(Graph::STATE_DEBLOOM_DONE))
        {
            if (graph._branchingKind != BRANCHING_NONE)
            {
                BranchingAlgorithm<span> branchingAlgo (
                    graph,
                    graph.getStorage(),
                    graph._branchingKind,
                    props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0,
                    props
                );
                executeAlgorithm (branchingAlgo, graph.getStorage(), props, graph._info);

                graph.setState(Graph::STATE_BRANCHING_DONE);

                /** We configure the variant. */
                data.setBranching (branchingAlgo.getBranchingCollection());
            }
        }

        /************************************************************/
        /*                    Post processing                       */
        /************************************************************/

        /** In case we choose another storage for the solid kmers, we have to update the graph state. */
        if (props->get(STR_URI_SOLID_KMERS) != 0)
        {
            data.setSolid (0);
            graph.unsetState (Graph::STATE_BANKCONVERTER_DONE);
            graph.unsetState (Graph::STATE_SORTING_COUNT_DONE);
        }

        /** We save library information in the root of the storage. */
        graph.getGroup().addProperty ("xml", string("\n") + LibraryInfo::getInfo().getXML());

        /** We save the state and kmer size at storage root level. */
        graph.getGroup().addProperty ("state",     Stringify::format("%d", graph._state));
        graph.getGroup().addProperty ("kmer_size", Stringify::format("%d", graph._kmerSize));

        /************************************************************/
        /*                        Clean up                          */
        /************************************************************/
    }

    /** Algorithm configuration. */
    void executeAlgorithm (Algorithm& algorithm, Storage& storage, IProperties* props, IProperties& info) const
    {
        algorithm.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));

        algorithm.run ();

        info.add (1, algorithm.getInfo());
        info.add (1, algorithm.getSystemInfo());

        /** We memorize information of the algorithm execution as a property of the corresponding group. */
        storage.getGroup(algorithm.getName()).addProperty("xml", string("\n") + algorithm.getInfo()->getXML());
    }
};


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
IOptionsParser* Graph::getOptionsParser (bool includeMandatory, bool enableMphf)
{

    /** We build the root options parser. */
    OptionsParser* parser = new OptionsParser ("graph");

    /** We add children parser to it (kmer count, bloom/debloom, branching). */
    parser->push_back (SortingCountAlgorithm<>::getOptionsParser(includeMandatory));
    parser->push_back (DebloomAlgorithm<>::getOptionsParser());
    parser->push_back (BranchingAlgorithm<>::getOptionsParser());

    /** We activate MPHF option only if available. */
    if (MPHF<char>::enabled)
    {
        IOptionsParser* parserEmphf  = new OptionsParser ("emphf");
        parserEmphf->push_back (new tools::misc::impl::OptionOneParam (STR_MPHF_TYPE, "mphf type ('none' or 'emphf')", false,  enableMphf ? "emphf":"none"));
        parser->push_back  (parserEmphf);
    }

    /** We create a "general options" parser. */
    IOptionsParser* parserGeneral  = new OptionsParser ("general");
    parserGeneral->push_front (new OptionOneParam (STR_INTEGER_PRECISION, "integers precision (0 for optimized value)", false, "0", false));
    parserGeneral->push_front (new OptionOneParam (STR_VERBOSE,           "verbosity level",      false, "1"  ));
    parserGeneral->push_front (new OptionOneParam (STR_NB_CORES,          "number of cores",      false, "0"  ));

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
Graph  Graph::create (bank::IBank* bank, const char* fmt, ...)
{
    IOptionsParser* parser = getOptionsParser (false);   LOCAL(parser);

    /** We build the command line from the format and the ellipsis. */
    std::string commandLine;
    char* buffer = 0;
    va_list args;
    va_start (args, fmt);
    int res = vasprintf (&buffer, fmt, args);
    va_end (args);
    if (buffer != NULL)  {  commandLine = buffer;  FREE (buffer);  }

    try
    {
        return  Graph (bank, parser->parseString(commandLine));
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
Graph  Graph::create (const char* fmt, ...)
{
    IOptionsParser* parser = getOptionsParser (true);   LOCAL (parser);

    /** We build the command line from the format and the ellipsis. */
    std::string commandLine;
    char* buffer = 0;
    va_list args;
    va_start (args, fmt);
    int res = vasprintf (&buffer, fmt, args);
    va_end (args);
    if (buffer != NULL)  {  commandLine = buffer;  FREE (buffer);  }

    try
    {
        return  Graph (parser->parseString(commandLine));
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
Graph::Graph (size_t kmerSize)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(kmerSize), _info("graph"),
      _state(Graph::STATE_INIT_DONE),
      _bloomKind(BLOOM_DEFAULT), _debloomKind(DEBLOOM_DEFAULT), _debloomImpl(DEBLOOM_IMPL_DEFAULT),
      _branchingKind(BRANCHING_STORED), _mphfKind(MPHF_NONE)
{
    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize);

    /** We configure the graph data from the storage content. */
    boost::apply_visitor (configure_visitor (*this, getStorage()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Graph (const std::string& uri)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph"), _name(System::file().getBaseName(uri))
{
    size_t precision = 0;

    /** We create a storage instance. */
    setStorage (StorageFactory(_storageMode).create (uri, false, false));

    /** We get some properties. */
    _state     = (Graph::StateMask) atol (getGroup().getProperty ("state").c_str());
    _kmerSize  =                    atol (getGroup().getProperty ("kmer_size").c_str());

    /** We get library information in the root of the storage. */
    string xmlString = getGroup().getProperty ("xml");
    stringstream ss; ss << xmlString;   IProperties* props = new Properties(); LOCAL(props);
    props->readXML (ss);  getInfo().add (1, props);

    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize);

    /** We configure the graph data from the storage content. */
    boost::apply_visitor (configure_visitor (*this, getStorage()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Graph (bank::IBank* bank, tools::misc::IProperties* params)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _state(Graph::STATE_INIT_DONE),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph")
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    size_t integerPrecision = params->getInt (STR_INTEGER_PRECISION);

    /** We get other user parameters. */
    parse (params->getStr(STR_BLOOM_TYPE),        _bloomKind);
    parse (params->getStr(STR_DEBLOOM_TYPE),      _debloomKind);
    parse (params->getStr(STR_DEBLOOM_IMPL),      _debloomImpl);
    parse (params->getStr(STR_BRANCHING_TYPE),    _branchingKind);

    /** This one is conditional. */
    if (params->get(STR_MPHF_TYPE)) {  parse (params->getStr(STR_MPHF_TYPE), _mphfKind); }
    else                            { _mphfKind = MPHF_NONE; }

    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize, integerPrecision);

    /** We build the graph according to the wanted precision. */
    boost::apply_visitor (build_visitor (*this, bank,params),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : this is mostly duplicated code with function above, TODO refactor?
*********************************************************************/
Graph::Graph (tools::misc::IProperties* params)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _state(Graph::STATE_INIT_DONE),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph")
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    size_t integerPrecision = params->getInt (STR_INTEGER_PRECISION);

    /** We get other user parameters. */
    parse (params->getStr(STR_BLOOM_TYPE),        _bloomKind);
    parse (params->getStr(STR_DEBLOOM_TYPE),      _debloomKind);
    parse (params->getStr(STR_DEBLOOM_IMPL),      _debloomImpl);
    parse (params->getStr(STR_BRANCHING_TYPE),    _branchingKind);

    /** This one is conditional. */
    if (params->get(STR_MPHF_TYPE)) {  parse (params->getStr(STR_MPHF_TYPE), _mphfKind); }
    else                            { _mphfKind = MPHF_NONE; }

    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize, integerPrecision);

    /** We build a Bank instance for the provided reads uri. */
    bank::IBank* bank = Bank::open (params->getStr(STR_URI_INPUT));

    /** We build the graph according to the wanted precision. */
    boost::apply_visitor (build_visitor (*this, bank,params),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Graph ()
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph"),
      _state(Graph::STATE_INIT_DONE),
      _bloomKind(BLOOM_DEFAULT),
      _debloomKind(DEBLOOM_DEFAULT), _debloomImpl(DEBLOOM_IMPL_DEFAULT), _branchingKind(BRANCHING_STORED), _mphfKind(MPHF_NONE)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Graph (const Graph& graph)
    : _storageMode(graph._storageMode), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(graph._kmerSize), _info("graph"), _name(graph._name)
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
Graph& Graph::operator= (const Graph& graph)
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
        _mphfKind        = graph._mphfKind;
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
Graph::~Graph ()
{
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
void Graph::remove ()
{
    getStorage().remove();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Graph::isBranching (const Node& node) const
{
    return (! (successors<Node>(node).size()==1 && predecessors<Node>(node).size() == 1));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool Graph::isSimple (const Edge& edge) const
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
bool Graph::isEdge (const Node& u, const Node& v) const
{
    bool result = false;

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
Node Graph::reverse (const Node& node) const
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
BranchingNode Graph::reverse (const BranchingNode& node) const
{
    BranchingNode result = node;
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
Edge Graph::reverse (const Edge& edge) const
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
size_t Graph::indegree  (const Node& node) const  {  return neighbors<Node> (node, DIR_INCOMING).size();   }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t Graph::outdegree (const Node& node) const  {  return neighbors<Node> (node, DIR_OUTCOMING).size();  }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t Graph::degree (const Node& node, Direction dir) const  {  return neighbors<Node> (node, dir).size();  }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Item, typename Functor>
struct getItems_visitor : public boost::static_visitor<Graph::Vector<Item> >    {

    const Node& source;  Direction direction;  Functor fct;

    getItems_visitor (const Node& aSource, Direction aDirection, Functor aFct) : source(aSource), direction(aDirection), fct(aFct) {}

    template<size_t span>  Graph::Vector<Item> operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        Graph::Vector<Item> items;

        size_t idx = 0;

        /** We get the specific typed value from the generic typed value. */
        const Type& sourceVal = source.kmer.get<Type>();

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        Type graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, kmerSize) );

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
                        fct (items, idx++, source.kmer, source.strand, Node::Value(forward), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        fct (items, idx++, source.kmer, source.strand, Node::Value(reverse), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
            }
        }

        if (direction & DIR_INCOMING)
        {
            /** IMPORTANT !!! Since we have hugely shift the nt value, we make sure to use a long enough integer. */
            for (u_int64_t nt=0; nt<4; nt++)
            {
                Type forward = ((graine >> 2 )  + ( Type(nt) << ((kmerSize-1)*2)) ) & mask; // previous kmer
                Type reverse = revcomp (forward, kmerSize);

                Nucleotide NT;

                if (source.strand == STRAND_FORWARD)  {  NT = (Nucleotide) (source.kmer[0]);  }
                else                                  {  NT = kmer::reverse ((Nucleotide) (source.kmer[(kmerSize-1)]));  }

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        fct (items, idx++, source.kmer, source.strand, Node::Value(forward), STRAND_FORWARD, NT, DIR_INCOMING);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        fct (items, idx++, source.kmer, source.strand, Node::Value(reverse), STRAND_REVCOMP, NT, DIR_INCOMING);
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

struct Functor_getEdges {  void operator() (
    Graph::Vector<Edge>& items,
    size_t               idx,
    const Node::Value&   kmer_from,
    kmer::Strand         strand_from,
    const Node::Value&   kmer_to,
    kmer::Strand         strand_to,
    kmer::Nucleotide     nt,
    Direction            dir
) const
{
    items[idx++].set (kmer_from, strand_from, kmer_to, strand_to, nt, dir);
}};

Graph::Vector<Edge> Graph::getEdges (const Node& source, Direction direction)  const
{

    return boost::apply_visitor (getItems_visitor<Edge,Functor_getEdges>(source, direction, Functor_getEdges()),  *(GraphDataVariant*)_variant);
}

/********************************************************************************/
struct Functor_getNodes {  void operator() (
    Graph::Vector<Node>&   items,
    size_t               idx,
    const Node::Value&   kmer_from,
    kmer::Strand         strand_from,
    const Node::Value&   kmer_to,
    kmer::Strand         strand_to,
    kmer::Nucleotide     nt,
    Direction            dir
) const
{
    items[idx++].set (kmer_to, strand_to);
}};

Graph::Vector<Node> Graph::getNodes (const Node& source, Direction direction)  const
{
    return boost::apply_visitor (getItems_visitor<Node,Functor_getNodes>(source, direction, Functor_getNodes()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Item, typename Functor>
struct getItemsCouple_visitor : public boost::static_visitor<Graph::Vector<pair<Item,Item> > >    {

    const Node& node1; const Node& node2;  Direction direction; Functor functor;

    getItemsCouple_visitor (const Node& node1, const Node& node2, Direction aDirection, Functor aFct)
        : node1(node1), node2(node2),  direction(aDirection), functor(aFct) {}

    template<size_t span>  Graph::Vector<pair<Item,Item> > operator() (const GraphData<span>& data) const
    {
        typedef typename Kmer<span>::Type  Type;

        size_t idx = 0;
        Graph::Vector < pair<Item,Item> > items;

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();

        /** We get the specific typed value from the generic typed value. */
        const Type& val1 = node1.kmer.get<Type>();
        const Type& val2 = node2.kmer.get<Type>();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
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
                        functor (p,
                            node1.kmer, node1.strand, Node::Value(forward1), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, Node::Value(forward2), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
                else if (isForwardMin1==true && isForwardMin2==false)
                {
                    if (data.contains (forward1) && data.contains (reverse2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        functor (p,
                            node1.kmer, node1.strand, Node::Value(forward1), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, Node::Value(reverse2), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
                else if (isForwardMin1==false && isForwardMin2==true)
                {
                    if (data.contains (reverse1) && data.contains (forward2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        functor (p,
                            node1.kmer, node1.strand, Node::Value(reverse1), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, Node::Value(forward2), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING
                        );
                    }
                }
                else if (isForwardMin1==false && isForwardMin2==false)
                {
                    if (data.contains (reverse1) && data.contains (reverse2))
                    {
                        pair<Item,Item>& p = items[idx++];
                        functor (p,
                            node1.kmer, node1.strand, Node::Value(reverse1), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING,
                            node2.kmer, node2.strand, Node::Value(reverse2), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING
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
struct Functor_getNodesCouple {  void operator() (
    pair<Node,Node>&   items,
    const Node::Value&   kmer_from1,
    kmer::Strand         strand_from1,
    const Node::Value&   kmer_to1,
    kmer::Strand         strand_to1,
    kmer::Nucleotide     nt1,
    Direction            dir1,
    const Node::Value&   kmer_from2,
    kmer::Strand         strand_from2,
    const Node::Value&   kmer_to2,
    kmer::Strand         strand_to2,
    kmer::Nucleotide     nt2,
    Direction            dir2
) const
{
    items.first.set  (kmer_to1, strand_to1);
    items.second.set (kmer_to2, strand_to2);
}};

Graph::Vector<std::pair<Node,Node> > Graph::getNodesCouple (const Node& node1, const Node& node2, Direction direction) const
{
    return boost::apply_visitor (getItemsCouple_visitor<Node,Functor_getNodesCouple>(node1, node2, direction, Functor_getNodesCouple()),  *(GraphDataVariant*)_variant);
}

/********************************************************************************/
struct Functor_getEdgesCouple {  void operator() (
        pair<Edge,Edge>&     items,
        const Node::Value&   kmer_from1,
        kmer::Strand         strand_from1,
        const Node::Value&   kmer_to1,
        kmer::Strand         strand_to1,
        kmer::Nucleotide     nt1,
        Direction            dir1,
        const Node::Value&   kmer_from2,
        kmer::Strand         strand_from2,
        const Node::Value&   kmer_to2,
        kmer::Strand         strand_to2,
        kmer::Nucleotide     nt2,
        Direction            dir2
) const
{
    items.first.set  (kmer_from1, strand_from1, kmer_to1, strand_to1, nt1, dir1);
    items.second.set (kmer_from2, strand_from2, kmer_to2, strand_to2, nt2, dir2);
}};

Graph::Vector<std::pair<Edge,Edge> > Graph::getEdgesCouple (const Node& node1, const Node& node2, Direction direction) const
{
    return boost::apply_visitor (getItemsCouple_visitor<Edge,Functor_getEdgesCouple>(node1, node2, direction, Functor_getEdgesCouple()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct buildNode_visitor : public boost::static_visitor<Node>    {

    const tools::misc::Data& data;  size_t offset;

    buildNode_visitor (const tools::misc::Data& aData, size_t aOffset) : data(aData), offset(aOffset)  {}

    template<size_t span>  Node operator() (const GraphData<span>& graphData) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::ModelCanonical::Kmer Kmer;

        Kmer kmer = graphData._model->getKmer (data, offset);
        return Node (Integer(kmer.value()), kmer.forward()==kmer.value() ? STRAND_FORWARD : STRAND_REVCOMP);
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
Node Graph::buildNode (const tools::misc::Data& data, size_t offset)  const
{
    return boost::apply_visitor (buildNode_visitor(data,offset),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Node Graph::buildNode (const char* sequence)  const
{
    Data data ((char*)sequence);

    return boost::apply_visitor (buildNode_visitor(data,0),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Vector<BranchingNode> Graph::getBranchingNodeNeighbors (const Node& source, Direction direction) const
{
    Graph::Vector<BranchingNode>  result;

    /** We get the neighbors of the source node. */
    Graph::Vector<Edge> neighbors = this->neighbors<Edge> (source, direction);

    /** We resize the result vector. */
    result.resize (neighbors.size());

    /** We loop over all the neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        /** We get a simple path iterator from the current neighbor. */
        Graph::Iterator<Edge> path = this->simplePath<Edge> (neighbors[i].to, direction);

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
Graph::Vector<BranchingEdge> Graph::getBranchingEdgeNeighbors (const Node& source, Direction direction) const
{
    Graph::Vector<BranchingEdge>  result;

    /** We get the neighbors of the source node. */
    Graph::Vector<Edge> neighbors = this->neighbors<Edge> (source, direction);

    /** We resize the result vector. */
    result.resize (neighbors.size());

    /** We loop over all the neighbors. */
    for (size_t i=0; i<neighbors.size(); i++)
    {
        DEBUG ((cout << "neighbor[" << i << "] " << this->toString(neighbors[i]) << endl));

        /** We get a simple path iterator from the current neighbor. */
        Graph::Iterator<Edge> path = this->simplePath<Edge> (neighbors[i].to, direction);

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
template<typename Item, typename Functor>
struct getItem_visitor : public boost::static_visitor<Item>    {

    const Node& source;  Direction direction;  Nucleotide nt; bool& exists; Functor fct;

    getItem_visitor (const Node& aSource, Direction aDirection, Nucleotide aNt, bool& aExists, Functor aFct)
        : source(aSource), direction(aDirection), nt(aNt), exists(aExists), fct(aFct) {}

    template<size_t span>  Item operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        Item item;

        /** We get the specific typed value from the generic typed value. */
        const Type& sourceVal = source.kmer.get<Type>();

        /** Shortcuts. */
        size_t      kmerSize = data._model->getKmerSize();
        const Type& mask     = data._model->getKmerMax();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        Type graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, kmerSize) );

        if (direction & DIR_OUTCOMING)
        {
            Type forward = ( (graine << 2 )  + nt) & mask;
            Type reverse = revcomp (forward, kmerSize);

            if (forward < reverse)
            {
                if (data.contains (forward))
                {
                    fct (item, source.kmer, source.strand, Node::Value(forward), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
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
                    fct (item, source.kmer, source.strand, Node::Value(reverse), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
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
            Type forward = ((graine >> 2 )  + ( Type(nt) << ((kmerSize-1)*2)) ) & mask; // previous kmer
            Type reverse = revcomp (forward, kmerSize);

            Nucleotide NT;

            if (source.strand == STRAND_FORWARD)  {  NT = (Nucleotide) (source.kmer[0]);  }
            else                                  {  NT = kmer::reverse ((Nucleotide) (source.kmer[(kmerSize-1)]));  }

            if (forward < reverse)
            {
                if (data.contains (forward))
                {
                    fct (item, source.kmer, source.strand, Node::Value(forward), STRAND_FORWARD, NT, DIR_INCOMING);
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
                    fct (item, source.kmer, source.strand, Node::Value(reverse), STRAND_REVCOMP, NT, DIR_INCOMING);
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

struct Functor_getNode {  void operator() (
    Node&                item,
    const Node::Value&   kmer_from,
    kmer::Strand         strand_from,
    const Node::Value&   kmer_to,
    kmer::Strand         strand_to,
    kmer::Nucleotide     nt,
    Direction            dir
) const
{
    item.set (kmer_to, strand_to);
}};


Node Graph::getNode (const Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{
    return boost::apply_visitor (getItem_visitor<Node,Functor_getNode>(source, dir, nt, exists, Functor_getNode()),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Vector<Edge> Graph::getEdgeValues (const Node::Value& kmer) const
{
    Node source (kmer);

    Graph::Vector<Edge> v1 = getEdges (source,          DIR_OUTCOMING);
    Graph::Vector<Edge> v2 = getEdges (reverse(source), DIR_OUTCOMING);
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
Graph::Vector<Node> Graph::getNodeValues (const Node::Value& kmer) const
{
    Node source (kmer);

    Graph::Vector<Node> v1 = getNodes (source,          DIR_OUTCOMING);
    Graph::Vector<Node> v2 = getNodes (reverse(source), DIR_OUTCOMING);

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
Graph::Vector<BranchingEdge> Graph::getBranchingEdgeValues (const Node::Value& kmer) const
{
    Node source (kmer);

    Graph::Vector<BranchingEdge> v1 = getBranchingEdgeNeighbors (source,          DIR_OUTCOMING);
    Graph::Vector<BranchingEdge> v2 = getBranchingEdgeNeighbors (reverse(source), DIR_OUTCOMING);
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
Graph::Vector<BranchingNode> Graph::getBranchingNodeValues (const Node::Value& kmer) const
{
    Node source (kmer);

    Graph::Vector<BranchingNode> v1 = getBranchingNodeNeighbors (source,          DIR_OUTCOMING);
    Graph::Vector<BranchingNode> v2 = getBranchingNodeNeighbors (reverse(source), DIR_OUTCOMING);

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
struct contains_visitor : public boost::static_visitor<bool>    {

    const Node& node;
    contains_visitor (const Node& aNode) : node(aNode) {}

    template<size_t span>  size_t operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        return data.contains (node.kmer.get<Type>());
    }
};

/********************************************************************************/
bool Graph::contains (const Node& item) const
{
    return boost::apply_visitor (contains_visitor(item),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename NodeType>  struct BranchingFilter
{
    const Graph& graph;
    BranchingFilter(const Graph& graph) : graph(graph) {}
    bool operator () (const NodeType& item) { return graph.isBranching(item); }
};

template<typename NodeType>
struct nodes_visitor : public boost::static_visitor<tools::dp::ISmartIterator<NodeType>*>
{
    const Graph& graph;
    nodes_visitor (const Graph& graph) : graph(graph) {}

    template<size_t span>  tools::dp::ISmartIterator<NodeType>* operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::ModelCanonical Model;
        typedef typename Kmer<span>::Type           Type;
        typedef typename Kmer<span>::Count          Count;

        class NodeIterator : public tools::dp::ISmartIterator<NodeType>
        {
        public:
            NodeIterator (tools::dp::Iterator<Count>* ref, u_int64_t nbItems)
                : _ref(0),  _rank(0), _isDone(true), _nbItems(nbItems)   {  setRef(ref);  this->_item->strand = STRAND_FORWARD; }

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
                    this->_rank ++;
                    this->_item->kmer      = _ref->item().value;
                    this->_item->abundance = _ref->item().abundance;
                }
            }

            /** \copydoc  Iterator::next */
            void next()
            {
                _ref->next();
                _isDone = _ref->isDone();
                if (!_isDone)
                {
                    this->_rank ++;
                    this->_item->kmer      = _ref->item().value;
                    this->_item->abundance = _ref->item().abundance;
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

        if (typeid(NodeType) == typeid(Node))
        {
            if (data._solid != 0)
            {
                return new NodeIterator (data._solid->iterator (), data._solid->getNbItems());
            }
            else
            {
                throw "Iteration impossible (no solid nodes available)";
            }
        }
        else if (typeid(NodeType) == typeid(BranchingNode))
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
                return new FilterIterator<NodeType,BranchingFilter<NodeType> > (
                    new NodeIterator (data._solid->iterator (), data._solid->getNbItems()),
                    BranchingFilter<NodeType> (graph)
                );
            }
            else
            {
                throw "Iteration impossible (no solid nor branching nodes available)";
            }
        }
        else {  throw "Invalid type";  }
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
Graph::Iterator<Node> Graph::getNodes () const
{
    return Graph::Iterator<Node> (boost::apply_visitor (nodes_visitor<Node>(*this),  *(GraphDataVariant*)_variant));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Iterator<BranchingNode> Graph::getBranchingNodes () const
{
    return Graph::Iterator<BranchingNode> (boost::apply_visitor (nodes_visitor<BranchingNode>(*this),  *(GraphDataVariant*)_variant));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct toString_node_visitor : public boost::static_visitor<std::string>    {

    const Node& node;
    toString_node_visitor (const Node& aNode) : node(aNode) {}

    template<size_t span>  std::string operator() (const GraphData<span>& data) const
    {
        /** Shortcut. */
        typedef typename Kmer<span>::Type Type;

        Type value = node.kmer.get<Type>();
        if (node.strand == STRAND_FORWARD)   {  return data._model->toString (value);  }
        else                                 {  return data._model->toString (data._model->reverse (value));  }
    }
};

std::string Graph::toString (const Node& node) const
{
    return boost::apply_visitor (toString_node_visitor(node),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
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

        Type value = node.kmer.get<Type>();

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

std::string Graph::debugString (const Node& node, kmer::Strand strand, int mode) const
{
    return boost::apply_visitor (debugString_node_visitor(node,strand,mode),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
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
               << debugString_node_visitor (edge.from, strand, 1) (data)
               << " ";
            if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
            else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }

            ss  << " "
                    << debugString_node_visitor (edge.to, strand, 1) (data)
               << "]";
        }
        else if (mode==2)
        {
            ss << "["
               << toString_node_visitor (edge.from) (data)
               << " ";
            if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
            else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }

            ss  << " "
                    << toString_node_visitor (edge.to) (data)
               << "]";
        }

        return ss.str();
    }
};

std::string Graph::debugString (const Edge& edge, kmer::Strand strand, int mode) const
{
    return boost::apply_visitor (debugString_edge_visitor(edge, strand, mode),  *(GraphDataVariant*)_variant);
}

std::string Graph::toString (const Edge& edge) const
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
std::string Graph::toString (const BranchingEdge& edge) const
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
int Graph::simplePathAvance (const Node& node, Direction dir, kmer::Nucleotide& nt) const
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
int Graph::simplePathAvance (const Node& node, Direction dir) const
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
int Graph::simplePathAvance (const Node& node, Direction dir, Edge& output) const
{
    Graph::Vector<Edge> neighbors = this->neighbors<Edge> (node, dir);

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

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<typename Item, typename Functor>
class AbstractSimplePathIterator : public tools::dp::ISmartIterator<Item>
{
public:

    AbstractSimplePathIterator (const Graph& graph, const Node& node, Direction dir, const Functor& update)
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
    const Graph&       _graph;
    Direction          _dir;
    u_int64_t          _rank;
    bool               _isDone;
    const Functor&     _update;
};

/** */
template<typename Functor>
class NodeSimplePathIterator : public AbstractSimplePathIterator<Node,Functor>
{
public:
    NodeSimplePathIterator (const Graph& graph, const Node& node, Direction dir, const Functor& update)
        : AbstractSimplePathIterator<Node,Functor> (graph, node, dir, update)  { *(this->_item) = node;  }
};

/** */
template<typename Functor>
class EdgeSimplePathIterator : public AbstractSimplePathIterator<Edge, Functor>
{
public:
    EdgeSimplePathIterator (const Graph& graph, const Node& node, Direction dir, const Functor& udpate)
        : AbstractSimplePathIterator<Edge, Functor> (graph, node, dir, udpate)
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
struct Functor_getSimpleNodeIterator {  void operator() (
    const Graph&         graph,
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

Graph::Iterator<Node> Graph::getSimpleNodeIterator (const Node& node, Direction dir) const
{
    return Graph::Iterator<Node> (new NodeSimplePathIterator <Functor_getSimpleNodeIterator> (*this, node, dir, Functor_getSimpleNodeIterator()));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct Functor_getSimpleEdgeIterator {  void operator() (
    const Graph&         graph,
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

Graph::Iterator<Edge> Graph::getSimpleEdgeIterator (const Node& node, Direction dir) const
{
    return Graph::Iterator<Edge> (new EdgeSimplePathIterator<Functor_getSimpleEdgeIterator>(*this, node, dir, Functor_getSimpleEdgeIterator()));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template <>
std::set<BranchingNode> Graph::neighbors (std::set<BranchingNode>::iterator first, std::set<BranchingNode>::iterator last) const
{
#if 0
    std::set<BranchingNode> result;
    for (auto it=first; it!=last; ++it)
    {
        Graph::Vector<BranchingNode> neighbors = this->neighbors<BranchingNode> (it->kmer);
        for (size_t i=0; i<neighbors.size(); i++)  { result.insert (neighbors[i]); }
    }
    return result;
#else

    static const size_t nbThread = 8;
    static const size_t nbThreshold = nbThread*1;

    std::set<BranchingNode> result;

    size_t nb = std::distance (first, last);

    if (nb >= nbThreshold)
    {
        std::set<BranchingNode>::iterator begin = first;
        std::set<BranchingNode>::iterator end;


        int nbPerThread = nb/nbThread;

        vector<pair<std::set<BranchingNode>::iterator,std::set<BranchingNode>::iterator> > iteratorPairs;

        class Cmd : public tools::dp::ICommand, public system::SmartPointer
        {
        public:
            Cmd (const Graph& graph, const pair<std::set<BranchingNode>::iterator,std::set<BranchingNode>::iterator>& range)
                : graph(graph), range(range)
            {
                result.reserve (std::distance(range.first,range.second)*8);
            }

            void execute ()
            {
                for (std::set<BranchingNode>::iterator it=range.first; it!=range.second; ++it)
                {
                    Graph::Vector<BranchingNode> neighbors = graph.neighbors<BranchingNode> (it->kmer);
                    for (size_t i=0; i<neighbors.size(); i++)  { result.push_back (neighbors[i]); }
                }
            }

            vector<BranchingNode>& get() { return result; }

        private:
            const Graph& graph;
            pair<std::set<BranchingNode>::iterator,std::set<BranchingNode>::iterator> range;
            vector<BranchingNode> result;
        };

        vector<tools::dp::ICommand*> cmds;

        while (end != last)
        {
            end = begin;
            advance (end, std::min (nbPerThread, (int)distance(end,last)));
            iteratorPairs.push_back (std::make_pair(begin, end) );

			tools::dp::ICommand* cmd = new Cmd (*this, std::make_pair(begin, end));
			cmd->use();
			cmds.push_back (cmd);
            begin = end;
        }

        tools::dp::impl::Dispatcher().dispatchCommands(cmds, 0);

        for (size_t i=0; i<cmds.size(); i++)
        {
            vector<BranchingNode>& current = ((Cmd*)cmds[i])->get();

            result.insert (current.begin(), current.end());
            cmds[i]->forget();
        }
    }
    else
    {
        for (std::set<BranchingNode>::iterator it=first; it!=last; ++it)
        {
            Graph::Vector<BranchingNode> neighbors = this->neighbors<BranchingNode> (it->kmer);
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
** REMARKS :
*********************************************************************/
struct mutate_visitor : public boost::static_visitor<Graph::Vector<Node> >    {

    const Node& node;  size_t idx;  int mode;

    mutate_visitor (const Node& node, size_t idx, int mode) : node(node), idx(idx), mode(mode){}

    template<size_t span>  Graph::Vector<Node> operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::Type           Type;
        typedef typename Kmer<span>::ModelCanonical Model;

        Graph::Vector<Node> result;
        size_t nbMutations = 0;

        size_t kmerSize = data._model->getKmerSize();

        Model* model = data._model;

        Type kmer = node.kmer.get<Type>();

        size_t idxReverse = kmerSize - 1 - idx;

        if (node.strand == STRAND_FORWARD)
        {
            Nucleotide ntInit = (Nucleotide) (node.kmer[idxReverse]);
            size_t     nt0    = mode==0 ? 0 : ntInit+1;

            Type resetMask = ((~(Type(3)<<(2*idxReverse))));

            for (size_t nt=nt0; nt<4; nt++)
            {
                Type direct = (kmer & resetMask) |   (Type(nt) << (2*idxReverse));
                Type rev    = model->reverse(direct);

                if (direct < rev)
                {
                    if (data.contains(direct))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = direct;
                        node.strand = STRAND_FORWARD;
                    }
                }
                else
                {
                    if (data.contains(rev))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = rev;
                        node.strand = STRAND_REVCOMP;
                    }
                }
            } /* end of or (size_t nt=nt0; */
        }
        else
        {
            Nucleotide ntInit = reverse ((Nucleotide) (node.kmer[idx]));
            size_t nt0 = mode==0 ? 0 : ntInit+1;

            Type resetMask = ((~(Type(3)<<(2*idx))));

            for (size_t nt=nt0; nt<4; nt++)
            {
                Type direct = (kmer & resetMask) +  (Type(reverse((Nucleotide)nt)) << (2*idx));
                Type rev    = model->reverse(direct);

                if (direct < rev)
                {
                    if (data.contains(direct))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = direct;
                        node.strand = STRAND_REVCOMP;
                    }
                }
                else
                {
                    if (data.contains(rev))
                    {
                        Node& node  = result[nbMutations++];
                        node.kmer   = rev;
                        node.strand = STRAND_FORWARD;
                    }
                }
            }
        }

        result.resize(nbMutations);

        return result;
    }
};

/** */
Graph::Vector<Node> Graph::mutate (const Node& node, size_t idx, int mode) const
{
    return boost::apply_visitor (mutate_visitor(node,idx,mode),  *(GraphDataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct getNT_visitor : public boost::static_visitor<Nucleotide>    {

    const Node& node;  size_t idx;

    getNT_visitor (const Node& node, size_t idx) : node(node), idx(idx) {}

    template<size_t span>  Nucleotide operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::ModelCanonical Model;
        size_t kmerSize = data._model->getKmerSize();

        if (node.strand == STRAND_FORWARD)  { return (Nucleotide) (node.kmer[kmerSize-1-idx]); }
        else                                { return reverse ((Nucleotide) (node.kmer[idx])); }
    }
};

/** */
Nucleotide Graph::getNT (const Node& node, size_t idx) const
{
    return boost::apply_visitor (getNT_visitor(node,idx),  *(GraphDataVariant*)_variant);
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/

// I don't understand fully this visitor pattern, was it needed for this method? -r
struct queryAbundance_visitor : public boost::static_visitor<int>    {

    const Node& node;

    queryAbundance_visitor (const Node& node) : node(node){}

    template<size_t span>  int operator() (const GraphData<span>& data) const
    {
        typedef typename Kmer<span>::Type  Type;
        unsigned char res = 0;

        /** We get the specific typed value from the generic typed value. */
        res = (*(data._abundance))[node.kmer.get<Type>()];

        return res;
    }
};

/** */
int Graph::queryAbundance (const Node& node) const
{
    return boost::apply_visitor (queryAbundance_visitor(node),  *(GraphDataVariant*)_variant);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
