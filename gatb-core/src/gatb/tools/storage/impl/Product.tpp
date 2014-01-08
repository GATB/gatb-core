/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/
namespace gatb  {  namespace core  {  namespace tools  {  namespace storage  {  namespace impl {
/********************************************************************************/

/********************************************************************************
     #####   #######  #        #        #######   #####   #######
    #     #  #     #  #        #        #        #     #     #
    #        #     #  #        #        #        #           #
    #        #     #  #        #        #####    #           #
    #        #     #  #        #        #        #           #
    #     #  #     #  #        #        #        #     #     #      ##
     #####   #######  #######  #######  #######   #####      #      ##
********************************************************************************/

/*********************************************************************
*********************************************************************/
template <class Item>
inline CollectionNode<Item>::CollectionNode (
    ProductFactory* factory, 
    ICell* parent, 
    const std::string& id, 
    collections::Collection<Item>* ref
)
    : _factory(factory), Cell(parent,id), collections::impl::CollectionAbstract<Item> (ref->bag(), ref->iterable()), _ref(0)
{
    /** We get a token on the referred Collection instance. */
    setRef(ref);
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline CollectionNode<Item>::~CollectionNode()
{
    /** We release the token of the Collection instance. */
    setRef(0);
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline void CollectionNode<Item>::remove ()
{
    /** We delegate the job to the referred Collection instance. */
    _ref->remove();
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline void CollectionNode<Item>::addProperty (const std::string& key, const std::string value)  
{  
    _ref->addProperty (key, value); 
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline std::string  CollectionNode<Item>::getProperty (const std::string& key) 
{  
    return _ref->getProperty (key); 
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline collections::Collection<Item>* CollectionNode<Item>::getRef ()  
{ 
    return _ref; 
}

/**********************************************************************
             #####   ######   #######  #     #  ######
            #     #  #     #  #     #  #     #  #     #
            #        #     #  #     #  #     #  #     #
            #  ####  ######   #     #  #     #  ######
            #     #  #   #    #     #  #     #  #
            #     #  #    #   #     #  #     #  #
             #####   #     #  #######   #####   #
**********************************************************************/

/*********************************************************************
*********************************************************************/
inline Group::Group (ProductFactory* factory, ICell* parent, const std::string& name) 
    : _factory(factory), Cell(parent, name) 
{
}

/*********************************************************************
*********************************************************************/
inline Group::~Group()
{
    /** We release the collections. */
    for (size_t i=0; i<_collections.size(); i++)  {  _collections[i]->forget();  }

    /** We release the partitions. */
    for (size_t i=0; i<_partitions.size(); i++)  { _partitions[i]->forget(); }

    /** We release the groups. */
    for (size_t i=0; i<_groups.size(); i++)  { _groups[i]->forget(); }
}

/*********************************************************************
*********************************************************************/
inline Group& Group::getGroup (const std::string& name)
{
    Group* group=0;  for (size_t i=0; !group && i<_groups.size(); i++)  {  if (_groups[i]->getId() == name)  { group = _groups[i]; }  }

    if (group == 0)
    {
        group = _factory->createGroup (this, name);
        _groups.push_back (group);
        group->use ();
    }
    return *group;
}

/*********************************************************************
*********************************************************************/
template <class Type>  
inline Partition<Type>&  Group::getPartition (const std::string& name, size_t nb)
{
    Partition<Type>* result = _factory->createPartition<Type> (this, name, nb);
    _partitions.push_back(result);
    result->use();
    return *result;
}

/*********************************************************************
*********************************************************************/
template <class Type>  
inline CollectionNode<Type>&  Group::getCollection (const std::string& name)
{
    CollectionNode<Type>* result = _factory->createCollection<Type> (this, name, 0);
    _collections.push_back (result);
    result->use ();
    return *result;
}

/*********************************************************************
*********************************************************************/
inline void Group::remove ()
{
    /** We remove the collections. */
    for (size_t i=0; i<_collections.size(); i++)  {  _collections[i]->remove ();  }

    /** We remove the partitions. */
    for (size_t i=0; i<_partitions.size(); i++)   { _partitions[i]->remove(); }

    /** We remove the children groups. */
    for (size_t i=0; i<_groups.size(); i++)       { _groups[i]->remove(); }
}

/**********************************************************************
######      #     ######   #######  ###  #######  ###  #######  #     #
#     #    # #    #     #     #      #      #      #   #     #  ##    #
#     #   #   #   #     #     #      #      #      #   #     #  # #   #
######   #     #  ######      #      #      #      #   #     #  #  #  #
#        #######  #   #       #      #      #      #   #     #  #   # #
#        #     #  #    #      #      #      #      #   #     #  #    ##
#        #     #  #     #     #     ###     #     ###  #######  #     #
**********************************************************************/

/*********************************************************************
*********************************************************************/
template<typename Type>
inline Partition<Type>::Partition (ProductFactory* factory, ICell* parent, const std::string& id, size_t nbCollections)
    : Group (factory, parent, id), _factory(factory), _typedCollections(nbCollections), _synchro(0)
{
    /** We create a synchronizer to be shared by the collections. */
    _synchro = system::impl::System::thread().newSynchronizer();

    /** We want to instantiate the wanted number of collections. */
    for (size_t i=0; i<_typedCollections.size(); i++)
    {
        /** We define the name of the current partition as a mere number. */
        std::stringstream ss;  ss << i;

        CollectionNode<Type>* result = _factory->createCollection<Type> (this, ss.str(), _synchro);

        /** We add the collection node to the dedicated vector and take a token for it. */
        (_typedCollections [i] = result)->use ();
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline Partition<Type>::~Partition ()
{
    /** We release the token for each collection node. */
    for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->forget ();  }

    delete _synchro;
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline size_t Partition<Type>::size()  const  
{ 
    return _typedCollections.size(); 
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline collections::Collection<Type>& Partition<Type>::operator[] (size_t idx)  
{  
    return  * _typedCollections[idx]->getRef();  
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void Partition<Type>::flush ()
{
    for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->flush();  }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void Partition<Type>::remove ()
{
    /** We remove each collection of this partition. */
    for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->remove ();  }

    /** We call the remove of the parent class. */
    Group::remove ();
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline PartitionCache<Type>::PartitionCache (Partition<Type>& ref, size_t nbItemsCache, system::ISynchronizer* synchro)
    :  _ref(ref), _nbItemsCache(nbItemsCache), _synchro(synchro), _cachedCollections(ref.size())
{
    /** We create the partition files. */
    for (size_t i=0; i<_cachedCollections.size(); i++)
    {
        _cachedCollections[i] = new collections::impl::CollectionCache<Type> (ref[i], nbItemsCache, synchro);
        _cachedCollections[i]->use ();
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline PartitionCache<Type>::PartitionCache (const PartitionCache<Type>& p)
    : _ref(p._ref), _nbItemsCache(p._nbItemsCache), _synchro(p._synchro), _cachedCollections(p.size())
{
    /** We create the partition files. */
    for (size_t i=0; i<_cachedCollections.size(); i++)
    {
        PartitionCache<Type>& pp = (PartitionCache<Type>&)p;
        _cachedCollections[i] = new collections::impl::CollectionCache<Type> (pp[i].getRef(), p._nbItemsCache, p._synchro);
        _cachedCollections[i]->use ();
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline PartitionCache<Type>::~PartitionCache ()
{
    flush ();
    for (size_t i=0; i<_cachedCollections.size(); i++)  {  _cachedCollections[i]->forget ();  }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline size_t  PartitionCache<Type>::size() const 
{ 
    return _cachedCollections.size(); 
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline collections::impl::CollectionCache<Type>&   PartitionCache<Type>::operator[] (size_t idx)  
{  
    return * _cachedCollections[idx];  
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void PartitionCache<Type>::flush ()   
{  
    for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->flush();  }  
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void PartitionCache<Type>::remove ()  
{  
    for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->remove ();  } 
}

/********************************************************************************
        ######   ######   #######  ######   #     #   #####   #######
        #     #  #     #  #     #  #     #  #     #  #     #     #
        #     #  #     #  #     #  #     #  #     #  #           #
        ######   ######   #     #  #     #  #     #  #           #
        #        #   #    #     #  #     #  #     #  #           #
        #        #    #   #     #  #     #  #     #  #     #     #
        #        #     #  #######  ######    #####    #####      #
********************************************************************************/

/*********************************************************************
*********************************************************************/
inline Product::Product (ProductMode_e mode, const std::string& name, bool autoRemove)
    : Cell(0, ""), _factory(0), _root(0), _autoRemove(autoRemove)  
{
    setFactory (new ProductFactory (mode));
}

inline Product::~Product ()
{
    setRoot    (0);
    setFactory (0);
}

/*********************************************************************
*********************************************************************/
inline Group* Product::getRoot ()
{
    if (_root == 0)  { setRoot    (_factory->createGroup (this, "")); }
    return _root;
}

/*********************************************************************
*********************************************************************/
inline Group& Product::operator() (const std::string name)
{
    if (name.empty())  { return *getRoot(); }
    else               { return getRoot ()->getGroup (name);  }
}

/*********************************************************************
*********************************************************************/
inline void Product::remove ()  
{  
    getRoot()->remove(); 
}

/*********************************************************************
*********************************************************************/
inline void Product::setFactory (ProductFactory* factory)  
{ 
    SP_SETATTR(factory); 
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/


#include <gatb/tools/storage/impl/ProductHDF5.hpp>
#include <gatb/tools/storage/impl/ProductFile.hpp>

/********************************************************************************/
namespace gatb  {  namespace core  {  namespace tools  {  namespace storage  {  namespace impl {
/********************************************************************************/

/********************************************************************************
        #######     #      #####   #######  #######  ######   #     #
        #          # #    #     #     #     #     #  #     #   #   #
        #         #   #   #           #     #     #  #     #    # #
        #####    #     #  #           #     #     #  ######      #
        #        #######  #           #     #     #  #   #       #
        #        #     #  #     #     #     #     #  #    #      #
        #        #     #   #####      #     #######  #     #     #
********************************************************************************/

/*********************************************************************
*********************************************************************/
inline Product* ProductFactory::createProduct (const std::string& name, bool deleteIfExist, bool autoRemove)
{
    switch (_mode)
    {
        case PRODUCT_HDF5:  return ProductHDF5Factory::createProduct (name, deleteIfExist, autoRemove);
        case PRODUCT_FILE:  return ProductFileFactory::createProduct (name, deleteIfExist, autoRemove);
        default:            throw system::Exception ("Unknown mode in ProductFactory::createProduct");
    }
}

/*********************************************************************
*********************************************************************/
inline Group* ProductFactory::createGroup (ICell* parent, const std::string& name)
{
    switch (_mode)
    {
        case PRODUCT_HDF5:  return ProductHDF5Factory::createGroup (parent, name);
        case PRODUCT_FILE:  return ProductFileFactory::createGroup (parent, name);
        default:            throw system::Exception ("Unknown mode in ProductFactory::createGroup");
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline Partition<Type>* ProductFactory::createPartition (ICell* parent, const std::string& name, size_t nb)
{
    switch (_mode)
    {
        case PRODUCT_HDF5:  return ProductHDF5Factory::createPartition<Type> (parent, name, nb);
        case PRODUCT_FILE:  return ProductFileFactory::createPartition<Type> (parent, name, nb);
        default:            throw system::Exception ("Unknown mode in ProductFactory::createPartition");
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline CollectionNode<Type>* ProductFactory::createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
{
    switch (_mode)
    {
        case PRODUCT_HDF5:  return ProductHDF5Factory::createCollection<Type> (parent, name, synchro);
        case PRODUCT_FILE:  return ProductFileFactory::createCollection<Type> (parent, name, synchro);
        default:            throw system::Exception ("Unknown mode in ProductFactory::createCollection");
    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

