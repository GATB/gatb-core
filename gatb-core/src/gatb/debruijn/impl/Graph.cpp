/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>
#include <gatb/debruijn/impl/BranchingAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/collections/impl/ContainerSet.hpp>

#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/DebloomAlgorithm.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

using namespace std;

using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

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
    typedef typename Kmer<span>::Model Model;
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

    /** Constructor. */
    GraphData () : _model(0), _solid(0), _bloom(0), _cfp(0), _branching(0) {}

    /** Destructor. */
    ~GraphData ()
    {
        setModel     (0);
        setSolid     (0);
        setBloom     (0);
        setCfp       (0);
        setBranching (0);
    }

    /** Constructor (copy). */
    GraphData (const GraphData& d) : _model(0), _solid(0), _bloom(0), _cfp(0), _branching(0)
    {
        setModel     (d._model);
        setSolid     (d._solid);
        setBloom     (d._bloom);
        setCfp       (d._cfp);
        setBranching (d._branching);
    }

    /** Assignment operator. */
    GraphData& operator= (const GraphData& d)
    {
        if (this != &d)
        {
            setModel     (d._model);
            setSolid     (d._solid);
            setBloom     (d._bloom);
            setCfp       (d._cfp);
            setBranching (d._branching);
        }
        return *this;
    }

    /** Required attributes. */
    Model*               _model;
    Collection<Count>*   _solid;
    Bloom<Type>*         _bloom;
    Container<Type>*     _cfp;
    Collection<Count>*   _branching;

    /** Setters. */
    void setModel       (Model*             model)      { SP_SETATTR (model);     }
    void setSolid       (Collection<Count>* solid)      { SP_SETATTR (solid);     }
    void setBloom       (Bloom<Type>*       bloom)      { SP_SETATTR (bloom);     }
    void setCfp         (Container<Type>*   cfp)        { SP_SETATTR (cfp);       }
    void setBranching   (Collection<Count>* branching)  { SP_SETATTR (branching); }

    /** Shortcut. */
    bool contains (const Type& item)  const  {  return _bloom->contains (item) && !_cfp->contains(item);  }
};

/********************************************************************************/

/* This definition is the basis for having a "generic" Graph class, ie. not relying on a template
 * parameter.
 *
 * This is done through a boost::variant; actually, we use a limited number of variant, corresponding
 * to some maximum kmer sizes.
 */
typedef boost::variant <
    GraphData<32>,
    GraphData<64>,
    GraphData<96>,
    GraphData<128>
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
static void setVariant (GraphDataVariant& data, size_t kmerSize)
{
    /** We convert the kmer size chosen by the user to the precision value. */
    size_t prec = 1 + kmerSize / 32;
    Integer::setType (prec);

    /** Here is the link between the kmer size (or precision) and the specific type
     * to be used for the variant. */
    switch (prec)
    {
        case 1:     data = GraphData<32> ();  break;
        case 2:     data = GraphData<64> ();  break;
        case 3:     data = GraphData<96> ();  break;
        case 4:     data = GraphData<128>();  break;
        default:    throw system::Exception ("Graph failure because of unhandled kmer precision %d", prec);
    }
}

/********************************************************************************/

/** This visitor is used to configure a GraphDataVariant object (ie configure its attributes).
 * The information source used to configure the variant is a kmer size and a storage.
 *
 * If the storage is null, only the kmer model is set.
 *
 * If the storage is not null, it is likely coming from a previous graph building (dsk, debloom, ...); we can
 * therefore configure the variant with the items coming from this storage.
 */
struct configure_visitor : public boost::static_visitor<>    {

    size_t   kmerSize;
    Storage* storage;

    configure_visitor (size_t aKmerSize, Storage* aStorage)  : kmerSize(aKmerSize), storage(aStorage) {}

    template<size_t span>  void operator() (GraphData<span>& data) const
    {
        /** We create the kmer model. */
        data.setModel (new typename Kmer<span>::Model (kmerSize));

        if (storage != NULL)
        {
            /** Shortcuts. */
            typedef typename Kmer<span>::Type  Type;
            typedef typename Kmer<span>::Count Count;

            /** We retrieve the needed information from the storage. */
            Collection<Count>*    solidIterable = & (*storage) ("dsk").    getCollection<Count >     ("solid");
            Iterable<Type>*       cFPKmers      = & (*storage) ("debloom").getCollection<Type>       ("cfp");
            Iterable<NativeInt8>* bloomArray    = & (*storage) ("debloom").getCollection<NativeInt8> ("bloom");
            Collection<Count>*    branching     = & (*storage) ("graph").  getCollection<Count >     ("branching");

            /** Some checks. */
            assert (solidIterable != 0);
            assert (cFPKmers      != 0);
            assert (bloomArray    != 0);
            assert (branching     != 0);

            /** We set the iterable for the solid kmers. */
            data.setSolid (solidIterable);

            /** We compute parameters for the Bloom filter. */
            double lg2 = log(2);
            float     NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
            size_t    nbHash         = (int)floorf (0.7*NBITS_PER_KMER);
            u_int64_t bloomSize      = (u_int64_t) (solidIterable->getNbItems() * NBITS_PER_KMER);

            /** We need a bloom builder for creating the Bloom filter. */
            kmer::impl::BloomBuilder<span> builder (bloomSize, nbHash);

            /** We check whether we have a stored Bloom filter. */
            if (bloomArray == 0)
            {
                /** We create Bloom filter and fill it with the solid kmers. */
                data.setBloom (builder.build (solidIterable->iterator()));
            }
            else
            {
                /** We load the Bloom filter from the specific dataset. */
                data.setBloom (builder.load (bloomArray));
            }

            /** We build the set of critical false positive kmers. */
            data.setCfp (new tools::collections::impl::ContainerSet<Type> (cFPKmers->iterator()));

            /** We set the branching container. */
            data.setBranching (branching);
        }
    }
};

/********************************************************************************/

/** This visitor is used to build a graph. In particular, the data variant of the graph will
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
 *  The data variant itself is configured after the debloom step. After that, the graph has enough
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

        size_t kmerSize = props->get(STR_KMER_SIZE)  ? props->getInt(STR_KMER_SIZE) : 27;
        size_t nks      = props->get(STR_NKS)        ? props->getInt(STR_NKS)       : 3;

        string output = props->get(STR_URI_OUTPUT) ?
            props->getStr(STR_URI_OUTPUT)   :
            (props->getStr(STR_URI_OUTPUT_DIR) + "/" + system::impl::System::file().getBaseName (bank->getId()));

        string binaryBankUri = System::file().getCurrentDirectory() + "/bank.bin";

        DEBUG ((cout << "builGraph for bank '" << bank->getId() << "'"
            << " kmerSize=" << kmerSize
            << " nks=" << nks
            << " output='" << output << "'"
            << endl
        ));

        /************************************************************/
        /*                       Storage creation                   */
        /************************************************************/
        graph.setStorage (StorageFactory(graph._storageMode).createStorage (output, true, false));

        /************************************************************/
        /*                         Bank conversion                  */
        /************************************************************/
        /** We create the binary bank. */
        BankConverterAlgorithm converter (bank, kmerSize, binaryBankUri);
        executeAlgorithm (converter, props, graph._info);

        /************************************************************/
        /*                         Sorting count                    */
        /************************************************************/
        /** We create a DSK instance and execute it. */
        SortingCountAlgorithm<span> sortingCount (
            graph._storage,
            converter.getResult(),
            kmerSize,
            nks,
            props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 1000,
            props->get(STR_MAX_DISK)   ? props->getInt(STR_MAX_DISK)   : 0,
            props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
        );
        executeAlgorithm (sortingCount, props, graph._info);

        /** We check that we got solid kmers. */
        if (sortingCount.getSolidKmers()->getNbItems() == 0)  {  throw "NO SOLID KMERS FOUND...";  }

        /************************************************************/
        /*                         Debloom                          */
        /************************************************************/
        /** We create a debloom instance and execute it. */
        DebloomAlgorithm<span> debloom (
            *graph._storage,
            sortingCount.getSolidKmers(),
            kmerSize,
            props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 1000,
            props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
        );
        executeAlgorithm (debloom, props, graph._info);

        /** We retrieve the bloom raw data. */
        Iterable<NativeInt8>* bloomArray = & graph.getStorage("debloom").getCollection<NativeInt8> ("bloom");

        /************************************************************/
        /*                 Variant configuration                    */
        /************************************************************/
        /** We configure the graph from the result of DSK and debloom. */
        boost::apply_visitor (configure_visitor (graph._kmerSize, graph._storage),  *(GraphDataVariant*)graph._variant);

        /************************************************************/
        /*                         Branching                        */
        /************************************************************/
        BranchingAlgorithm<span> branchingAlgo (graph, & graph.getStorage("graph").getCollection<Count> ("branching"));
        executeAlgorithm (branchingAlgo, props, graph._info);

        /************************************************************/
        /*                    Post processing                       */
        /************************************************************/
        /** We add metadata to some collections. */
        sortingCount.getSolidKmers()->addProperty ("properties", sortingCount.getInfo()->getXML());
        debloom.getCriticalKmers()  ->addProperty ("properties", debloom.     getInfo()->getXML());

        /** We add a special collection for global metadata (in particular the kmer size). */
        Collection<NativeInt8>* metadata = & graph.getStorage().getCollection<NativeInt8> ("metadata");
        NativeInt8 kmerSizeData[] = { kmerSize };
        metadata->insert (kmerSizeData, 1);

        metadata->addProperty ("properties", graph.getInfo().getXML());

        /************************************************************/
        /*                        Clean up                          */
        /************************************************************/
        /** We can get rid of the binary bank. */
        System::file().remove (binaryBankUri);
    }

    /** Algorithm configuration. */
    void executeAlgorithm (Algorithm& algorithm, IProperties* props, IProperties& info) const
    {
        algorithm.getInput()->add (0, STR_VERBOSE, props->getStr(STR_VERBOSE));

        algorithm.execute();

        info.add (1, algorithm.getInfo());
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
tools::misc::impl::OptionsParser Graph::getOptionsParser (bool includeMandatory)
{
    tools::misc::impl::OptionsParser parser;

    if (includeMandatory == true)
    {
        parser.push_back (new tools::misc::impl::OptionOneParam (STR_URI_INPUT, "reads file", true ));
    }

    parser.push_back (new tools::misc::impl::OptionOneParam (STR_KMER_SIZE,       "size of a kmer",                       false,  "31"    ));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_URI_OUTPUT,      "output file",                          false));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_URI_OUTPUT_DIR,  "output directory",                     false,  "."));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_VERBOSE,         "verbosity level",                      false,  "1"));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_MAX_MEMORY,      "max memory",                           false, "1000"));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_MAX_DISK,        "max disk",                             false, "0"));
    parser.push_back (new tools::misc::impl::OptionOneParam (STR_NB_CORES,        "nb cores (0 for all)",                 false, "0"));
    parser.push_back (new tools::misc::impl::OptionNoParam  (STR_HELP,            "help",                                 false));

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
    OptionsParser parser = getOptionsParser (false);

    /** We build the command line from the format and the ellipsis. */
    std::string commandLine;
    char* buffer = 0;
    va_list args;
    va_start (args, fmt);
    vasprintf (&buffer, fmt, args);
    va_end (args);
    if (buffer != NULL)  {  commandLine = buffer;  free (buffer);  }

    return  Graph (bank, parser.parse(commandLine));
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
    OptionsParser parser = getOptionsParser (true);

    /** We build the command line from the format and the ellipsis. */
    std::string commandLine;
    char* buffer = 0;
    va_list args;
    va_start (args, fmt);
    vasprintf (&buffer, fmt, args);
    va_end (args);
    if (buffer != NULL)  {  commandLine = buffer;  free (buffer);  }

    return  Graph (parser.parse(commandLine));
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
      _variant(new GraphDataVariant()), _kmerSize(kmerSize), _info("graph")
{
    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize);

    /** We configure the graph data from the storage content. */
    boost::apply_visitor (configure_visitor (_kmerSize, _storage),  *(GraphDataVariant*)_variant);
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
    setStorage (StorageFactory(_storageMode).createStorage (uri, false, false));

    /** We retrieve the type of kmers to be used from the storage. */
    Collection<NativeInt8>* metadata = & getStorage().getCollection<NativeInt8> ("metadata");
    gatb::core::tools::dp::Iterator<NativeInt8>* itData = metadata->iterator();  LOCAL (itData);
    itData->first(); if (!itData->isDone())  { _kmerSize = itData->item(); }

    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize);

    /** We configure the graph data from the storage content. */
    boost::apply_visitor (configure_visitor (_kmerSize, _storage),  *(GraphDataVariant*)_variant);
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
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph")
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize);

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
Graph::Graph (tools::misc::IProperties* params)
    : _storageMode(PRODUCT_MODE_DEFAULT), _storage(0),
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph")
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    /** We configure the data variant according to the provided kmer size. */
    setVariant (*((GraphDataVariant*)_variant), _kmerSize);

    /** We build a Bank instance for the provided reads uri. */
    bank::IBank* bank = new BankFasta (params->getStr(STR_URI_INPUT));

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
      _variant(new GraphDataVariant()), _kmerSize(0), _info("graph")
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
        _kmerSize = graph._kmerSize;

        _storageMode = graph._storageMode;
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
    DEBUG ((cout << "Graph::remove  NOT IMPLEMENTED..." << endl));
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
        const Type& mask     = data._model->getMask();

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
Graph::Vector<Edge> Graph::getEdges (const Node& source, Direction direction)  const
{
    struct Functor {  void operator() (
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

    return boost::apply_visitor (getItems_visitor<Edge,Functor>(source, direction, Functor()),  *(GraphDataVariant*)_variant);
}

/********************************************************************************/
Graph::Vector<Node> Graph::getNodes (const Node& source, Direction direction)  const
{
    struct Functor {  void operator() (
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

    return boost::apply_visitor (getItems_visitor<Node,Functor>(source, direction, Functor()),  *(GraphDataVariant*)_variant);
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
        typedef typename Kmer<span>::Type Type;

        pair<Type,bool> kmer = graphData._model->getKmer (data, offset);
        return Node (Integer(kmer.first), kmer.second ? STRAND_FORWARD : STRAND_REVCOMP);
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
        const Type& mask     = data._model->getMask();

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

Node Graph::getNode (const Node& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{
    struct Functor {  void operator() (
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

    return boost::apply_visitor (getItem_visitor<Node,Functor>(source, dir, nt, exists, Functor()),  *(GraphDataVariant*)_variant);
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
template<typename NodeType>
struct nodes_visitor : public boost::static_visitor<tools::dp::ISmartIterator<NodeType>*>    {

    template<size_t span>  tools::dp::ISmartIterator<NodeType>* operator() (const GraphData<span>& data) const
    {
        /** Shortcuts. */
        typedef typename Kmer<span>::Model Model;
        typedef typename Kmer<span>::Type  Type;
        typedef typename Kmer<span>::Count Count;

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

             if (typeid(NodeType) == typeid(Node))           {  return new NodeIterator (data._solid->iterator (),     data._solid->getNbItems());  }
        else if (typeid(NodeType) == typeid(BranchingNode))  {  return new NodeIterator (data._branching->iterator (), data._branching->getNbItems());  }
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
    return Graph::Iterator<Node> (boost::apply_visitor (nodes_visitor<Node>(),  *(GraphDataVariant*)_variant));
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
    return Graph::Iterator<BranchingNode> (boost::apply_visitor (nodes_visitor<BranchingNode>(),  *(GraphDataVariant*)_variant));
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
               << debugString_node_visitor (edge.from, strand, 0) (data)
               << " ";
            if (edge.direction == DIR_OUTCOMING)  {  ss <<  "--"  << ascii(edge.nt) << "-->";  }
            else                                  {  ss <<  "<--" << ascii(edge.nt) << "--";   }

            ss  << " "
                    << debugString_node_visitor (edge.to, strand, 0) (data)
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
Graph::Iterator<Node> Graph::getSimpleNodeIterator (const Node& node, Direction dir) const
{
    struct Functor {  void operator() (
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

    return Graph::Iterator<Node> (new NodeSimplePathIterator <Functor> (*this, node, dir, Functor()));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Graph::Iterator<Edge> Graph::getSimpleEdgeIterator (const Node& node, Direction dir) const
{
    struct Functor {  void operator() (
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

    return Graph::Iterator<Edge> (new EdgeSimplePathIterator<Functor>(*this, node, dir, Functor()));
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
                for (auto it=range.first; it!=range.second; ++it)
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
        for (auto it=first; it!=last; ++it)
        {
            Graph::Vector<BranchingNode> neighbors = this->neighbors<BranchingNode> (it->kmer);
            for (size_t i=0; i<neighbors.size(); i++)  { result.insert (neighbors[i]); }
        }
    }

    return result;
#endif
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
