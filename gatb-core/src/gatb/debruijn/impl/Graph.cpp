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


#define DEBUG(a)  //printf a

/********************************************************************************/

/** IMPORTANT !!!
 * We want to have the same behavior than the original minia, in particular the order
 * of the incoming neighbors. Since we don't use exactly the same way to compute
 * incoming neighbors, we have to restore the same order with a little reordering trick.
 */
u_int64_t incomingTable[] = { 2, 3, 0, 1 };
//u_int64_t incomingTable[] = { 0, 1, 2, 3 };
bool hack = true;


/********************************************************************************/
namespace gatb {  namespace core {  namespace debruijn {  namespace impl {
/********************************************************************************/

#if 0
#define ProductFactoryLocal tools::collections::impl::ProductFileFactory
#else
#define ProductFactoryLocal tools::collections::impl::ProductHDF5Factory
#endif

/********************************************************************************/

/** We define a structure that holds all the necessary stuff for implementing the graph API.
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

/** This definition is the basis for having a "generic" Graph class, ie. not relying on a template
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

/** This class has 3 main parts:
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
            system::impl::System::file().getBaseName (bank->getId());

        string binaryBankUri = System::file().getCurrentDirectory() + "/bank.bin";

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
    tools::misc::impl::OptionsParser parser;
    parser.add (new tools::misc::impl::OptionOneParam (STR_KMER_SIZE, "size of a kmer",                       true            ));
    parser.add (new tools::misc::impl::OptionOneParam (STR_NKS,       "abundance threshold for solid kmers",  false,  "3"     ));

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
    : _product(0), _variant(new DataVariant()), _kmerSize(0), _info("graph")
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
Graph::Graph (const Graph& graph) : _product(0), _variant(new DataVariant()), _kmerSize(graph._kmerSize), _info("graph")
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
    DEBUG (("Graph::remove  NOT IMPLEMENTED...\n"));
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
    begin = node;    for (std::vector<Node> nodes ; (nodes = predecessors<Node> (begin)).size() == 1;  begin = nodes[0])  {}
    end   = node;    for (std::vector<Node> nodes ; (nodes = successors  <Node> (end  )).size() == 1;  end   = nodes[0])  {}
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
Node Graph::getNode (const tools::misc::Data& data, size_t offset)  const
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

    Nucleotide NT;

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
struct getItems_visitor : public boost::static_visitor<std::vector<Item> >    {

    const Node& source;  Direction direction;  Functor fct;

    getItems_visitor (const Node& source, Direction direction, Functor fct) : source(source), direction(direction), fct(fct) {}

    template<typename T>  std::vector<Item> operator() (const Data<T>& data) const
    {
        std::vector<Item> items(8);

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
            for (u_int64_t k=0; k<ARRAY_SIZE(incomingTable); k++)
            {
                u_int64_t nt = incomingTable[k];

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
std::vector<Edge> Graph::getEdges (const Node& source, Direction direction)  const
{
    struct Functor {  void operator() (
        std::vector<Edge>& items,
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


    std::vector<Edge> result =  boost::apply_visitor (getItems_visitor<Edge,Functor>(source, direction, Functor()),  *(DataVariant*)_variant);

#ifdef GRAPH_CHECK
    for (size_t i=0; i<result.size(); i++)
    {
        Edge& edge = result[i];

        string from = toString (edge.from, STRAND_FORWARD, 1);
        string to   = toString (edge.to,   STRAND_FORWARD, 1);
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
std::vector<Node> Graph::getNodes (const Node& source, Direction direction)  const
{
    struct Functor {  void operator() (
        std::vector<Node>&   items,
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
std::vector<Edge> Graph::getEdgeValues (const Node::Value& kmer) const
{
    Node source (kmer);

    std::vector<Edge> v1 = getEdges (source,          DIR_OUTCOMING);
    std::vector<Edge> v2 = getEdges (reverse(source), DIR_OUTCOMING);

    /** We reverse the edges of v2. */
    for (size_t i=0; i<v2.size(); i++)
    {
        swap (v2[i].from, v2[i].to);
        v2[i].direction = impl::reverse (v2[i].direction);
    }

    v1.insert (v1.end(), v2.begin(), v2.end());

    return v1;
}

/********************************************************************************/
std::vector<Node> Graph::getNodeValues (const Node::Value& kmer) const
{
    Node source (kmer);

    std::vector<Node> v1 = getNodes (source,          DIR_OUTCOMING);
    std::vector<Node> v2 = getNodes (reverse(source), DIR_OUTCOMING);

    v1.insert (v1.end(), v2.begin(), v2.end());

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
                : _ref(0),  _isDone(true), _nbItems(nbItems)   {  setRef(ref);  this->_item->strand = STRAND_FORWARD; }

            ~NodeIterator ()  { setRef(0);   }

            /** \copydoc  Iterator::first */
            void first()
            {
                _ref->first();
                _isDone = _ref->isDone();

                if (!_isDone)
                {
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
            u_int64_t getNbItems () const { return _nbItems; }

        private:
            tools::dp::Iterator<kmer::Kmer<T> >* _ref;
            void setRef (tools::dp::Iterator<kmer::Kmer<T> >* ref)  { SP_SETATTR(ref); }

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

        if (mode==0)    {  ss << "[ " << kmerStr <<  " / " << strandStr << "]";  }
        else            {  ss << kmerStr; }

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

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
