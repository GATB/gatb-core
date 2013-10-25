/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/debruijn/impl/Graph.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/collections/impl/ContainerSet.hpp>

#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

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

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

#define DEBUG(a)  //printf a

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
    Data () : _model(0), _solid(0), _bloom(0), _cfp(0) {}

    /** Destructor. */
    ~Data ()
    {
        setModel (0);
        setSolid (0);
        setBloom (0);
        setCfp   (0);
    }

    /** Constructor (copy). */
    Data (const Data& d) : _model(0), _solid(0), _bloom(0), _cfp(0)
    {
        setModel (d._model);
        setSolid (d._solid);
        setBloom (d._bloom);
        setCfp   (d._cfp);
    }

    /** Assignment operator. */
    Data& operator= (const Data& d)
    {
        if (this != &d)
        {
            setModel (d._model);
            setSolid (d._solid);
            setBloom (d._bloom);
            setCfp   (d._cfp);
        }
        return *this;
    }

    /** Required attributes. */
    Model<T>*                   _model;
    Collection<kmer::Kmer<T> >* _solid;
    Bloom<T>*                   _bloom;
    Container<T>*               _cfp;

    /** Setters. */
    void setModel (Model<T>* model)                     { SP_SETATTR (model); }
    void setSolid (Collection<kmer::Kmer<T> >* solid)   { SP_SETATTR (solid); }
    void setBloom (Bloom<T>* bloom)                     { SP_SETATTR (bloom); }
    void setCfp   (Container<T>* cfp)                   { SP_SETATTR (cfp);   }

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
        configureVariant (graph, kmerSize, sortingCount.getSolidKmers(), debloom.getCriticalKmers(), bloomArray);

        /** We can get rid of the binary bank. */
        System::file().remove (binaryBankUri);

        /** We add metadata to some collections. */
        sortingCount.getSolidKmers()->addProperty ("properties", sortingCount.getInfo()->getXML());
        debloom.getCriticalKmers()  ->addProperty ("properties", debloom.     getInfo()->getXML());

        /** We add a special collection for global metadata (in particular the kmer size). */
        Collection<NativeInt8>* metadata = & graph.getProduct().getCollection<NativeInt8> ("metadata");
        NativeInt8 kmerSizeData[] = { kmerSize };
        metadata->insert (kmerSizeData, 1);
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
            & (*product) ("debloom").getCollection<NativeInt8>     ("bloom")
        );
    }

    /********************************************************************************/
    static void configureVariant (
        Graph& graph,
        size_t kmerSize,
        tools::collections::Collection<Kmer<T> >* solidIterable,
        tools::collections::Iterable<T>*          cFPKmers,
        tools::collections::Iterable<NativeInt8>* bloomArray
    )
    {
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
Graph::Graph (bank::IBank* bank, tools::misc::IProperties* params)
    : _product(0), _variant(new DataVariant())
{
    /** We get the kmer size from the user parameters. */
    size_t kmerSize = params->getInt (STR_KMER_SIZE);

    /** We get the Integer precision. */
    int precision = 1 + kmerSize / 32;
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
    : _product(0), _variant(new DataVariant())
{
    /** We get the kmer size from the user parameters. */
    size_t kmerSize = params->getInt (STR_KMER_SIZE);

    /** We build a Bank instance for the provided reads uri. */
    bank::IBank* bank = new Bank (params->getStr(STR_URI_INPUT));

    /** We get the Integer precision. */
    int precision = 1 + kmerSize / 32;
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
    : _product(0), _variant(new DataVariant())
{
    size_t kmerSize  = 0;
    size_t precision = 0;

    /** We create a product instance. */
    setProduct (ProductFactoryLocal::createProduct (uri, false, false));

    /** We retrieve the type of kmers to be used from the product. */
    Collection<NativeInt8>* metadata = & getProduct().getCollection<NativeInt8> ("metadata");
    Iterator<NativeInt8>* itData = metadata->iterator();  LOCAL (itData);
    itData->first(); if (!itData->isDone())  { kmerSize = itData->item(); }

    /** We set the precision of Integer to be used. */
    precision = 1 + kmerSize / 32;
    Integer::setType (precision);

    /** We load the graph. */
    switch (precision)
    {
        case 1: GraphFactoryImpl<LargeInt<1> >::loadGraph (*this, uri, kmerSize);  break;
        case 2: GraphFactoryImpl<LargeInt<2> >::loadGraph (*this, uri, kmerSize);  break;
        case 3: GraphFactoryImpl<LargeInt<3> >::loadGraph (*this, uri, kmerSize);  break;
        case 4: GraphFactoryImpl<LargeInt<4> >::loadGraph (*this, uri, kmerSize);  break;
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
Graph::~Graph ()
{
    /** We release resources. */
    setProduct (0);
    delete (DataVariant*)_variant;
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
    NodeSet nodes (*this);

    begin = node;    for (size_t nbNodes=0 ; (nbNodes = getPredecessors (begin, nodes)) == 1;  begin = nodes[0])  {}
    end   = node;    for (size_t nbNodes=0 ; (nbNodes = getSuccessors   (end,   nodes)) == 1;  end   = nodes[0])  {}
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct getEdges_visitor : public boost::static_visitor<size_t>    {

    const Node& source;  EdgeSet& edges;  Direction direction;
    getEdges_visitor (const Node& source, EdgeSet& edges, Direction direction) : source(source), edges(edges), direction(direction) {}

    template<typename T>  size_t operator() (const Data<T>& data) const
    {
        size_t idx = 0;

        /** We get the specific typed value from the generic typed value. */
        const T& sourceVal = source.kmer.value.get<T>();

        /** Shortcuts. */
        size_t   span = data._model->getSpan();
        const T& mask = data._model->getMask();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        T graine = ((source.strand == STRAND_FORWARD) ?  sourceVal :  revcomp (sourceVal, span) );

        if (direction == DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ( (graine << 2 )  + nt) & mask;
                T reverse = revcomp (forward, span);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        edges[idx++].set (source.getGraph(), source.kmer.value, source.strand, Type(forward), STRAND_FORWARD, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        edges[idx++].set (source.getGraph(), source.kmer.value, source.strand, Type(reverse), STRAND_REVCOMP, (Nucleotide)nt, DIR_OUTCOMING);
                    }
                }
            }
        }
        else if (direction == DIR_INCOMING)
        {
            /** IMPORTANT !!! Since we have hugely shift the nt value, we make sure to use a long enough integer. */
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ((graine >> 2 )  + ( nt << ((span-1)*2)) ) & mask; // previous kmer
                T reverse = revcomp (forward, span);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        edges[idx++].set (source.getGraph(), Type(forward), STRAND_FORWARD, source.kmer.value, source.strand, (Nucleotide)nt, DIR_INCOMING);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        edges[idx++].set (source.getGraph(), Type(reverse), STRAND_REVCOMP, source.kmer.value, source.strand, (Nucleotide)nt, DIR_INCOMING);
                    }
                }
            }
        }

        else  {  throw system::Exception ("Unknown direction for getting edges");  }

        /** We update the size of the container according to the number of found items. */
        edges.setSize (idx);

        /** We return the number of found items. */
        return idx;
    }
};

/********************************************************************************/
size_t Graph::getEdges (const Node& source, EdgeSet& edges, Direction direction)  const
{
    return boost::apply_visitor (getEdges_visitor(source, edges, direction),  *(DataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct getNodes_visitor : public boost::static_visitor<size_t>    {

    const Node& source;  NodeSet& nodes;  Direction direction;
    getNodes_visitor (const Node& source, NodeSet& nodes, Direction direction) : source(source), nodes(nodes), direction(direction) {}

    template<typename T>  size_t operator() (const Data<T>& data) const
    {
        size_t idx = 0;

        /** We get the specific typed value from the generic typed value. */
        const T& sourceVal = source.kmer.value.get<T>();

        /** Shortcuts. */
        size_t   span = data._model->getSpan();
        const T& mask = data._model->getMask();

        // the kmer we're extending may be actually a revcomp sequence in the bidirected debruijn graph node
        T graine = ((source.strand == STRAND_FORWARD) ?  sourceVal : revcomp (sourceVal, span));

        if (direction == DIR_OUTCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ( (graine << 2 )  + nt) & mask;    // next kmer
                T reverse = revcomp (forward, span);

                if (forward < reverse)
                {
                    if (data.contains (forward))
                    {
                        nodes[idx++].set (source.getGraph(), Type(forward), STRAND_FORWARD);
                    }
                }
                else
                {
                    if (data.contains (reverse))
                    {
                        nodes[idx++].set (source.getGraph(), Type(reverse), STRAND_REVCOMP);
                    }
                }
            }
        }
        else if (direction == DIR_INCOMING)
        {
            for (u_int64_t nt=0; nt<4; nt++)
            {
                T forward = ((graine >> 2 )  + ( nt <<  ((span-1)*2)) ) & mask; // previous kmer
                T reverse = revcomp (forward, span);

                if (forward < reverse)  {  if (data.contains (forward))  {  nodes[idx++].set (source.getGraph(), Type(forward), STRAND_FORWARD);  }  }
                else                    {  if (data.contains (reverse))  {  nodes[idx++].set (source.getGraph(), Type(reverse), STRAND_REVCOMP);  }  }
            }
        }

        else  {  throw system::Exception ("Unknown direction for getting nodes");  }

        /** We update the size of the container according to the number of found items. */
        nodes.setSize (idx);

        /** We return the number of found items. */
        return idx;
    }
};

/********************************************************************************/
size_t Graph::getNodes (const Node& source, NodeSet& nodes, Direction direction)  const
{
    return boost::apply_visitor (getNodes_visitor(source, nodes, direction),  *(DataVariant*)_variant);
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
        return data.contains (node.kmer.value.get<T>());
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
struct nodes_visitor : public boost::static_visitor<INodeIterator*>    {

    template<typename T>  INodeIterator* operator() (const Data<T>& data) const
    {
        class NodeIterator : public INodeIterator
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
                    this->_item->kmer.value     = _ref->item().value;
                    this->_item->kmer.abundance = _ref->item().abundance;
                }
            }

            /** \copydoc  Iterator::next */
            void next()
            {
                _ref->next();
                _isDone = _ref->isDone();
                if (!_isDone)
                {
                    this->_item->kmer.value     = _ref->item().value;
                    this->_item->kmer.abundance = _ref->item().abundance;
                }
            }

            /** \copydoc  Iterator::isDone */
            bool isDone() { return _isDone;  }

            /** \copydoc  Iterator::item */
            Node& item ()  {  return *(this->_item);  }

            /** */
            void setItem (Node& i)
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

        return new NodeIterator (data._solid->iterator (), data._solid->getNbItems());
    }
};


/********************************************************************************/
INodeIterator* Graph::nodes () const
{
    DEBUG (("Graph::nodes \n"));
    return boost::apply_visitor (nodes_visitor(),  *(DataVariant*)_variant);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
struct toString_visitor : public boost::static_visitor<std::string>    {

    const Node& node;  kmer::Strand& strand; int& mode;
    toString_visitor (const Node& node, kmer::Strand strand, int mode) : node(node), strand(strand), mode(mode) {}

    template<typename T>  std::string operator() (const Data<T>& data) const
    {
        std::stringstream ss;

        /** We set the strings for the kmer and the strand. */
        std::string kmerStr;
        std::string strandStr;

        T value = node.kmer.value.get<T>();

        if (strand == STRAND_ALL || node.strand == strand)
        {
            kmerStr = data._model->toString (value);
            strandStr = (node.strand==STRAND_FORWARD ? "FORWARD" : "REVCOMP");
        }
        else
        {
            T reverse = data._model->reverse (value);
            kmerStr = data._model->toString (reverse);
            strandStr = (node.strand==STRAND_FORWARD ? "REVCOMP" : "FORWARD");
        }

        if (mode==0)    {  ss << "[ " << kmerStr <<  "  strand=" << strandStr << "  abund=" << node.kmer.abundance << "]";  }
        else            {  ss << kmerStr; }

        return ss.str();

    }
};

std::string Graph::toString (const Node& node, kmer::Strand strand, int mode) const
{
    return boost::apply_visitor (toString_visitor(node,strand,mode),  *(DataVariant*)_variant);
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
