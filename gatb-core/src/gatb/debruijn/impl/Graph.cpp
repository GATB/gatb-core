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

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

//using namespace gatb::core::tools::dp;
//using namespace gatb::core::tools::dp::impl;

#undef NDEBUG
#include <cassert>

#define DEBUG(a)  //a

/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

#if 0
#define ProductFactoryLocal tools::collections::impl::ProductFileFactory
#else
#define ProductFactoryLocal tools::collections::impl::ProductHDF5Factory
#endif

/********************************************************************************/

/* We define a structure that holds all the necessary stuff for implementing the graph API.
 *  Here, the structure is templated by the required Integer class.
 *
 *  This structure is the basis for defining a boost::variant with all required integer
 *  template instantiations.
 */
template<typename T>
struct Data
{
    /** Constructor. */
    Data () : _model(0), _solid(0), _bloom(0), _cfp(0), _branching(0) {}

    /** Destructor. */
    ~Data ()
    {
        setModel     (0);
        setSolid     (0);
        setBloom     (0);
        setCfp       (0);
        setBranching (0);
    }

    /** Constructor (copy). */
    Data (const Data& d) : _model(0), _solid(0), _bloom(0), _cfp(0), _branching(0)
    {
        setModel     (d._model);
        setSolid     (d._solid);
        setBloom     (d._bloom);
        setCfp       (d._cfp);
        setBranching (d._branching);
    }

    /** Assignment operator. */
    Data& operator= (const Data& d)
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
    Model<T>*                   _model;
    Collection<kmer::Kmer<T> >* _solid;
    Bloom<T>*                   _bloom;
    Container<T>*               _cfp;
    Collection<kmer::Kmer<T> >* _branching;

    /** Setters. */
    void setModel       (Model<T>* model)                       { SP_SETATTR (model);     }
    void setSolid       (Collection<kmer::Kmer<T> >* solid)     { SP_SETATTR (solid);     }
    void setBloom       (Bloom<T>* bloom)                       { SP_SETATTR (bloom);     }
    void setCfp         (Container<T>* cfp)                     { SP_SETATTR (cfp);       }
    void setBranching   (Collection<kmer::Kmer<T> >* branching) { SP_SETATTR (branching); }

    /** Shortcut. */
    bool contains (const T& item)  const  {  return _bloom->contains (item) && !_cfp->contains(item);  }
};

/********************************************************************************/

/* This definition is the basis for having a "generic" Graph class, ie. not relying on a template
 * parameter.
 *
 * This is done through a boost::variant; actually, we use a limited number of variant, corresponding
 * to some maximum kmer size.
 */
typedef boost::variant <
    Data<tools::math::LargeInt<1> >,
    Data<tools::math::LargeInt<2> >,
    Data<tools::math::LargeInt<3> >,
    Data<tools::math::LargeInt<4> >
>  DataVariant;

/********************************************************************************/

/* This class has 3 main parts:
 *
 *      1) buildGraph: Graph creation from user parameters and save in filesystem
 *
 *      2) loadGraph:  Graph load from the filesystem for a given uri
 *
 *      3) configureVariant: set the variant attribute of the Graph class according
 *         to the real T template type.
 */
template<typename T>
class GraphFactoryImpl
{
public:

    /********************************************************************************/
    static void buildGraph (Graph& graph, bank::IBank* bank, tools::misc::IProperties* props)
    {
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
        /*                       Product creation                   */
        /************************************************************/
        graph.setProduct (ProductFactoryLocal::createProduct (output, true, false));

        /************************************************************/
        /*                         Bank conversion                  */
        /************************************************************/
        /** We create the binary bank. */
        BankConverterAlgorithm converter (bank, kmerSize, binaryBankUri);
        graph.executeAlgorithm (converter, props, graph._info);

        /************************************************************/
        /*                         Sorting count                    */
        /************************************************************/
        /** We create a DSK instance and execute it. */
        SortingCountAlgorithm<ProductFactoryLocal, T> sortingCount (
            graph._product,
            converter.getResult(),
            kmerSize,
            nks,
            props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 1000,
            props->get(STR_MAX_DISK)   ? props->getInt(STR_MAX_DISK)   : 0,
            props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
        );
        graph.executeAlgorithm (sortingCount, props, graph._info);

        /** We check that we got solid kmers. */
        if (sortingCount.getSolidKmers()->getNbItems() == 0)  {  throw "NO SOLID KMERS FOUND...";  }

        /************************************************************/
        /*                         Debloom                          */
        /************************************************************/
        /** We create a debloom instance and execute it. */
        DebloomAlgorithm<ProductFactoryLocal, T> debloom (
            *graph._product,
            sortingCount.getSolidKmers(),
            kmerSize,
            props->get(STR_MAX_MEMORY) ? props->getInt(STR_MAX_MEMORY) : 1000,
            props->get(STR_NB_CORES)   ? props->getInt(STR_NB_CORES)   : 0
        );
        graph.executeAlgorithm (debloom, props, graph._info);

        /** We retrieve the bloom raw data. */
        Iterable<NativeInt8>* bloomArray = & graph.getProduct("debloom").getCollection<NativeInt8> ("bloom");

        /** We configure the graph from the result of DSK and debloom. */
        configureVariant (graph, kmerSize, sortingCount.getSolidKmers(), debloom.getCriticalKmers(), bloomArray, 0);

        /** We can get rid of the binary bank. */
        System::file().remove (binaryBankUri);

        /** We add metadata to some collections. */
        sortingCount.getSolidKmers()->addProperty ("properties", sortingCount.getInfo()->getXML());
        debloom.getCriticalKmers()  ->addProperty ("properties", debloom.     getInfo()->getXML());

        /** We add a special collection for global metadata (in particular the kmer size). */
        Collection<NativeInt8>* metadata = & graph.getProduct().getCollection<NativeInt8> ("metadata");
        NativeInt8 kmerSizeData[] = { kmerSize };
        metadata->insert (kmerSizeData, 1);

        /************************************************************/
        /*                         Branching                        */
        /************************************************************/
        BranchingAlgorithm<ProductFactoryLocal, T> branchingAlgo (graph, & graph.getProduct("graph").getCollection<kmer::Kmer<T> > ("branching"));
        graph.executeAlgorithm (branchingAlgo, props, graph._info);
    }

    /********************************************************************************/
    static void loadGraph (Graph& graph, const std::string uri, size_t kmerSize)
    {
        /** We create a product instance. */
        Product<ProductFactoryLocal>* product = graph._product;

        /** We configure the graph with the aggregated information. */
        configureVariant (
            graph,
            kmerSize,
            & (*product) ("dsk").    getCollection<kmer::Kmer<T> > ("solid"),
            & (*product) ("debloom").getCollection<T >             ("cfp"),
            & (*product) ("debloom").getCollection<NativeInt8>     ("bloom"),
            & (*product) ("graph").  getCollection<kmer::Kmer<T> > ("branching")
        );
    }

    /********************************************************************************/
    static void configureVariant (
        Graph& graph,
        size_t kmerSize,
        tools::collections::Collection<Kmer<T> >* solidIterable,
        tools::collections::Iterable<T>*          cFPKmers,
        tools::collections::Iterable<NativeInt8>* bloomArray,
        tools::collections::Collection<Kmer<T> >* branching
    )
    {
        /** We set the kmer size. */
        graph._kmerSize = kmerSize;

        /** We create a templated data and configure it with instances needed for
         * implementing the graph API. */
        Data<T> data;

        /** We create the kmer model. */
        data.setModel (new kmer::impl::Model<T> (kmerSize));

        /** We set the iterable for the solid kmers. */
        data.setSolid (solidIterable);

        /** We compute parameters for the Bloom filter. */
        double lg2 = log(2);
        float     NBITS_PER_KMER = log (16*kmerSize*(lg2*lg2))/(lg2*lg2);
        size_t    nbHash         = (int)floorf (0.7*NBITS_PER_KMER);
        u_int64_t bloomSize      = (u_int64_t) (solidIterable->getNbItems() * NBITS_PER_KMER);

        /** We need a bloom builder for creating the Bloom filter. */
        kmer::impl::BloomBuilder<T> builder (bloomSize, nbHash);

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
        data.setCfp (new tools::collections::impl::ContainerSet<T> (cFPKmers->iterator()));

        /** We set the branching container. */
        data.setBranching (branching);

        /** THE BIG IDEA IS HERE... We set the Graph variant with the specific T data we just configured.
         *  From this moment, the Graph instance has the needed resources for providing the services,
         *  without any explicit reference to a LargeInt<T> type. Actually, the "lost" T type information
         *  will be retrieved through the use of boost::static_visitor class.
         */
        *((DataVariant*)graph._variant) = data;
    }

    /********************************************************************************/
    static void configureVariant (Graph& graph, size_t kmerSize)
    {
        /** We create a templated data and configure it with instances needed for
         * implementing the graph API. */
        Data<T> data;

        /** We create the kmer model. */
        data.setModel (new kmer::impl::Model<T> (kmerSize));

        /** THE BIG IDEA IS HERE... We set the Graph variant with the specific T data we just configured.
         *  From this moment, the Graph instance has the needed resources for providing the services,
         *  without any explicit reference to a LargeInt<T> type. Actually, the "lost" T type information
         *  will be retrieved through the use of boost::static_visitor class.
         */
        *((DataVariant*)graph._variant) = data;
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
        parser.add (new tools::misc::impl::OptionOneParam (STR_URI_INPUT, "reads file", true ));
    }

    parser.add (new tools::misc::impl::OptionOneParam (STR_KMER_SIZE,       "size of a kmer",                       false,  "31"    ));
    parser.add (new tools::misc::impl::OptionOneParam (STR_NKS,             "abundance threshold for solid kmers",  false,  "3"     ));
    parser.add (new tools::misc::impl::OptionOneParam (STR_URI_OUTPUT,      "output file",                          false));
    parser.add (new tools::misc::impl::OptionOneParam (STR_URI_OUTPUT_DIR,  "output directory",                     false,  "."));
    parser.add (new tools::misc::impl::OptionNoParam  (STR_VERBOSE,         "verbose",                              false));
    parser.add (new tools::misc::impl::OptionOneParam (STR_MAX_MEMORY,      "max memory",                           false, "1000"));
    parser.add (new tools::misc::impl::OptionOneParam (STR_MAX_DISK,        "max disk",                             false, "0"));
    parser.add (new tools::misc::impl::OptionOneParam (STR_NB_CORES,        "nb cores (0 for all)",                 false, "0"));
    parser.add (new tools::misc::impl::OptionNoParam  (STR_HELP,            "help",                                 false));

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
: _product(0), _variant(new DataVariant()), _kmerSize(kmerSize), _info("graph")
{
    /** We get the Integer precision. */
    int precision = 1 + _kmerSize / 32;
    Integer::setType (precision);

    /** We build the graph according to the wanted precision. */
    switch (precision)
    {
        case 1: GraphFactoryImpl<LargeInt<1> >::configureVariant (*this, kmerSize);  break;
        case 2: GraphFactoryImpl<LargeInt<2> >::configureVariant (*this, kmerSize);  break;
        case 3: GraphFactoryImpl<LargeInt<3> >::configureVariant (*this, kmerSize);  break;
        case 4: GraphFactoryImpl<LargeInt<4> >::configureVariant (*this, kmerSize);  break;
        default:   throw system::Exception ("Graph failure because of unhandled kmer precision %d", precision);
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
Graph::Graph (bank::IBank* bank, tools::misc::IProperties* params)
    : _product(0), _variant(new DataVariant()), _kmerSize(0), _info("graph")
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    /** We get the Integer precision. */
    int precision = 1 + _kmerSize / 32;
    Integer::setType (precision);

    /** We build the graph according to the wanted precision. */
    switch (precision)
    {
        case 1: GraphFactoryImpl<LargeInt<1> >::buildGraph (*this, bank, params);  break;
        case 2: GraphFactoryImpl<LargeInt<2> >::buildGraph (*this, bank, params);  break;
        case 3: GraphFactoryImpl<LargeInt<3> >::buildGraph (*this, bank, params);  break;
        case 4: GraphFactoryImpl<LargeInt<4> >::buildGraph (*this, bank, params);  break;
        default:   throw system::Exception ("Graph failure because of unhandled kmer precision %d", precision);
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
Graph::Graph (tools::misc::IProperties* params)
    : _product(0), _variant(new DataVariant()), _kmerSize(0), _info("graph")
{
    /** We get the kmer size from the user parameters. */
    _kmerSize = params->getInt (STR_KMER_SIZE);

    /** We build a Bank instance for the provided reads uri. */
    bank::IBank* bank = new Bank (params->getStr(STR_URI_INPUT));

    /** We get the Integer precision. */
    int precision = 1 + _kmerSize / 32;
    Integer::setType (precision);

    /** We build the graph according to the wanted precision. */
    switch (precision)
    {
        case 1: GraphFactoryImpl<LargeInt<1> >::buildGraph (*this, bank, params);  break;
        case 2: GraphFactoryImpl<LargeInt<2> >::buildGraph (*this, bank, params);  break;
        case 3: GraphFactoryImpl<LargeInt<3> >::buildGraph (*this, bank, params);  break;
        case 4: GraphFactoryImpl<LargeInt<4> >::buildGraph (*this, bank, params);  break;
        default:   throw system::Exception ("Graph failure because of unhandled kmer precision %d", precision);
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
Graph::Graph (const std::string& uri)
    : _product(0), _variant(new DataVariant()), _kmerSize(0), _info("graph"), _name(System::file().getBaseName(uri))
{
    size_t precision = 0;

    /** We create a product instance. */
    setProduct (ProductFactoryLocal::createProduct (uri, false, false));

    /** We retrieve the type of kmers to be used from the product. */
    Collection<NativeInt8>* metadata = & getProduct().getCollection<NativeInt8> ("metadata");
    gatb::core::tools::dp::Iterator<NativeInt8>* itData = metadata->iterator();  LOCAL (itData);
    itData->first(); if (!itData->isDone())  { _kmerSize = itData->item(); }

    /** We set the precision of Integer to be used. */
    precision = 1 + _kmerSize / 32;
    Integer::setType (precision);

    /** We load the graph. */
    switch (precision)
    {
        case 1: GraphFactoryImpl<LargeInt<1> >::loadGraph (*this, uri, _kmerSize);  break;
        case 2: GraphFactoryImpl<LargeInt<2> >::loadGraph (*this, uri, _kmerSize);  break;
        case 3: GraphFactoryImpl<LargeInt<3> >::loadGraph (*this, uri, _kmerSize);  break;
        case 4: GraphFactoryImpl<LargeInt<4> >::loadGraph (*this, uri, _kmerSize);  break;
        default:   throw system::Exception ("Graph failure because of unhandled kmer precision %d", precision);
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
Graph::Graph () : _product(0), _variant(new DataVariant()), _kmerSize(0), _info("graph")
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
Graph::Graph (const Graph& graph) : _product(0), _variant(new DataVariant()), _kmerSize(graph._kmerSize), _info("graph"), _name(graph._name)
{
    setProduct (graph._product);

    if (graph._variant)  {  *((DataVariant*)_variant) = *((DataVariant*)graph._variant);  }
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

        setProduct (graph._product);

        if (graph._variant)  {  *((DataVariant*)_variant) = *((DataVariant*)graph._variant);  }
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
    setProduct (0);
    if (_variant)  {  delete (DataVariant*)_variant;  }
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
void Graph::getNearestBranchingRange (const Node& node, Node& begin, Node& end) const
{
    begin = node;    for (Graph::Vector<Node> nodes ; (nodes = predecessors<Node> (begin)).size() == 1;  begin = nodes[0])  {}
    end   = node;    for (Graph::Vector<Node> nodes ; (nodes = successors  <Node> (end  )).size() == 1;  end   = nodes[0])  {}
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
Node Graph::buildNode (const tools::misc::Data& data, size_t offset)  const
{
    /** We create a kmer model. */
    Model<Integer> model (this->getKmerSize());

    pair<Integer,bool> kmer = model.getKmer (data, offset);

    Strand strand = kmer.second ? STRAND_FORWARD : STRAND_REVCOMP;

    return Node (kmer.first, strand);
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

    getItems_visitor (const Node& source, Direction direction, Functor fct) : source(source), direction(direction), fct(fct) {}

    template<typename T>  Graph::Vector<Item> operator() (const Data<T>& data) const
    {
        Graph::Vector<Item> items;

        size_t idx = 0;

        /** We get the specific typed value from the generic typed value. */
        const T& sourceVal = source.kmer.get<T>();

        /** Shortcuts. */
        size_t   span = data._model->getSpan();
        const T& mask = data._model->getMask();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        T graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, span) );

        if (direction & DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ( (graine << 2 )  + nt) & mask;
                T reverse = revcomp (forward, span);

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
                T forward = ((graine >> 2 )  + ( T(nt) << ((span-1)*2)) ) & mask; // previous kmer
                T reverse = revcomp (forward, span);

                Nucleotide NT;

                if (source.strand == STRAND_FORWARD)  {  NT = (Nucleotide) (source.kmer[0]);  }
                else                                  {  NT = kmer::reverse ((Nucleotide) (source.kmer[(span-1)]));  }

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


    Graph::Vector<Edge> result =  boost::apply_visitor (getItems_visitor<Edge,Functor>(source, direction, Functor()),  *(DataVariant*)_variant);

#ifdef GRAPH_CHECK
    for (size_t i=0; i<result.size(); i++)
    {
        Edge& edge = result[i];

        string from = toString (edge.from);
        string to   = toString (edge.to);
        string nt; nt   = (ascii(edge.nt));
        string check;
        if (edge.direction==DIR_OUTCOMING)
        {
            check = from + nt;  check = check.substr (1, check.size()-1);
            assert (check == to);
        }
        else
        {
            check = to   + nt;  check = check.substr (1, check.size()-1);
            assert (check == from);
        }
    }
#endif

    return result;
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

    return boost::apply_visitor (getItems_visitor<Node,Functor>(source, direction, Functor()),  *(DataVariant*)_variant);
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

    getItem_visitor (const Node& source, Direction direction, Nucleotide nt, bool& exists, Functor fct)
        : source(source), direction(direction), nt(nt), exists(exists), fct(fct) {}

    template<typename T>  Item operator() (const Data<T>& data) const
    {
        Item item;

        /** We get the specific typed value from the generic typed value. */
        const T& sourceVal = source.kmer.get<T>();

        /** Shortcuts. */
        size_t   span = data._model->getSpan();
        const T& mask = data._model->getMask();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        T graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, span) );

        if (direction & DIR_OUTCOMING)
        {
            T forward = ( (graine << 2 )  + nt) & mask;
            T reverse = revcomp (forward, span);

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
            T forward = ((graine >> 2 )  + ( T(nt) << ((span-1)*2)) ) & mask; // previous kmer
            T reverse = revcomp (forward, span);

            Nucleotide NT;

            if (source.strand == STRAND_FORWARD)  {  NT = (Nucleotide) (source.kmer[0]);  }
            else                                  {  NT = kmer::reverse ((Nucleotide) (source.kmer[(span-1)]));  }

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

    return boost::apply_visitor (getItem_visitor<Node,Functor>(source, dir, nt, exists, Functor()),  *(DataVariant*)_variant);
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
    contains_visitor (const Node& node) : node(node) {}

    template<typename T>  size_t operator() (const Data<T>& data) const
    {
        return data.contains (node.kmer.get<T>());
    }
};

/********************************************************************************/
bool Graph::contains (const Node& item) const
{
    return boost::apply_visitor (contains_visitor(item),  *(DataVariant*)_variant);
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

    template<typename T>  tools::dp::ISmartIterator<NodeType>* operator() (const Data<T>& data) const
    {
        class NodeIterator : public tools::dp::ISmartIterator<NodeType>
        {
        public:
            NodeIterator (tools::dp::Iterator<kmer::Kmer<T> >* ref, u_int64_t nbItems)
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
            tools::dp::Iterator<kmer::Kmer<T> >* _ref;
            void setRef (tools::dp::Iterator<kmer::Kmer<T> >* ref)  { SP_SETATTR(ref); }

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
    return Graph::Iterator<Node> (boost::apply_visitor (nodes_visitor<Node>(),  *(DataVariant*)_variant));
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
    return Graph::Iterator<BranchingNode> (boost::apply_visitor (nodes_visitor<BranchingNode>(),  *(DataVariant*)_variant));
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
    toString_node_visitor (const Node& node) : node(node) {}

    template<typename T>  std::string operator() (const Data<T>& data) const
    {
        T value = node.kmer.get<T>();
        if (node.strand == STRAND_FORWARD)   {  return data._model->toString (value);  }
        else                                 {  return data._model->toString (data._model->reverse (value));  }
    }
};

std::string Graph::toString (const Node& node) const
{
    return boost::apply_visitor (toString_node_visitor(node),  *(DataVariant*)_variant);
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

    const Node& node;  kmer::Strand& strand; int& mode;
    debugString_node_visitor (const Node& node, kmer::Strand strand, int mode) : node(node), strand(strand), mode(mode) {}

    template<typename T>  std::string operator() (const Data<T>& data) const
    {
        std::stringstream ss;

        /** We set the strings for the kmer and the strand. */
        std::string kmerStr;
        std::string strandStr;

        T value = node.kmer.get<T>();

        if (strand == STRAND_ALL || node.strand == strand)
        {
            kmerStr = data._model->toString (value);
            strandStr = (node.strand==STRAND_FORWARD ? "FWD" : "REV");
        }
        else
        {
            T reverse = data._model->reverse (value);
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
    return boost::apply_visitor (debugString_node_visitor(node,strand,mode),  *(DataVariant*)_variant);
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
    debugString_edge_visitor (const Edge& edge, Strand strand, int mode=0) : edge(edge), strand(strand), mode(mode) {}

    template<typename T>  std::string operator() (const Data<T>& data) const
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
    return boost::apply_visitor (debugString_edge_visitor(edge, strand, mode),  *(DataVariant*)_variant);
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
void Graph::executeAlgorithm (Algorithm& algorithm, IProperties* props, IProperties& info)
{
    string bargraph = props->get(STR_VERBOSE) ?  "2" : "0";

    algorithm.getInput()->add (0, STR_PROGRESS_BAR, bargraph);

    algorithm.execute();

    info.add (1, algorithm.getInfo());
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
