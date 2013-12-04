/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Product.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/impl/Cell.hpp>

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/CollectionAbstract.hpp>
#include <gatb/tools/collections/impl/CollectionCache.hpp>
#include <gatb/tools/collections/impl/CollectionFile.hpp>

#include <gatb/tools/misc/api/IProperty.hpp>

#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <map>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** Forward declarations. */
template <typename Factory>  class Group;
template <typename Factory, typename Type>  class Partition;
template <typename Type>  class CollectionNode;

/********************************************************************************
     #####   #######  #        #        #######   #####   #######
    #     #  #     #  #        #        #        #     #     #
    #        #     #  #        #        #        #           #
    #        #     #  #        #        #####    #           #
    #        #     #  #        #        #        #           #
    #     #  #     #  #        #        #        #     #     #      ##
     #####   #######  #######  #######  #######   #####      #      ##
********************************************************************************/
/** \brief Class that add some features to the Collection interface
 *
 * The CollectionNode has two aspects:
 *      - this is a Collection
 *      - this is a Node
 *
 * The idea was not to have a Collection interface extending the INode interface
 * => we introduce the CollectionNode for this.
 */
template <class Item>
class CollectionNode : public dp::impl::Cell, public CollectionAbstract<Item>
{
public:

    /** Constructor.
     *  The idea is to use a referred collection for:
     *      - its bag part
     *      - its iterable part
     *      - its remove part
     * \param[in] parent : parent node
     * \param[in] id  : identifier of the collection to be created
     * \param[in] ref : referred collection.
     */
    CollectionNode (dp::ICell* parent, const std::string& id, Collection<Item>* ref)
        : dp::impl::Cell(parent,id), CollectionAbstract<Item> (ref->bag(), ref->iterable()), _ref(0)
    {
        /** We get a token on the referred Collection instance. */
        setRef(ref);
    }

    /** Destructor. */
    virtual ~CollectionNode()
    {
        /** We release the token of the Collection instance. */
        setRef(0);
    }

    /** \copydoc dp::ICell::remove */
    void remove ()
    {
        /** We delegate the job to the referred Collection instance. */
        _ref->remove();
    }

    /** */
    void addProperty (const std::string& key, const std::string value)  {  _ref->addProperty (key, value); }

    /** */
    std::string getProperty (const std::string& key) {  return _ref->getProperty (key); }

    /**  */
    Collection<Item>* getRef ()  { return _ref; }

private:

    Collection<Item>* _ref;
    void setRef (Collection<Item>* ref)  { SP_SETATTR(ref); }
};

/**********************************************************************
             #####   ######   #######  #     #  ######
            #     #  #     #  #     #  #     #  #     #
            #        #     #  #     #  #     #  #     #
            #  ####  ######   #     #  #     #  ######
            #     #  #   #    #     #  #     #  #
            #     #  #    #   #     #  #     #  #
             #####   #     #  #######   #####   #
**********************************************************************/

/** */
template <typename Factory>
class Group : public dp::impl::Cell
{
public:

    /** Constructor. */
    Group (dp::ICell* parent, const std::string& name) : dp::impl::Cell(parent, name) {}

    /** Destructor. */
    ~Group()
    {
        /** We release the collections. */
        for (size_t i=0; i<_collections.size(); i++)  {  _collections[i]->forget();  }

        /** We release the partitions. */
        for (size_t i=0; i<_partitions.size(); i++)  { _partitions[i]->forget(); }

        /** We release the groups. */
        for (size_t i=0; i<_groups.size(); i++)  { _groups[i]->forget(); }
    }

    /** */
    Group<Factory>& getGroup (const std::string& name)
    {
        Group* group=0;  for (size_t i=0; !group && i<_groups.size(); i++)  {  if (_groups[i]->getId() == name)  { group = _groups[i]; }  }

        if (group == 0)
        {
            group = Factory::createGroup (this, name);
            _groups.push_back (group);
            group->use ();
        }
        return *group;
    }

    /** */
    template <class Type>  Partition<Factory,Type>& getPartition (const std::string& name, size_t nb)
    {
        /** See http://stackoverflow.com/questions/15572415/expected-primary-expression-before-in-g-but-not-in-microsoft-compiler */
        Partition<Factory,Type>* result = Factory::template createPartition<Type> (this, name, nb);
        _partitions.push_back(result);
        result->use();
        return *result;
    }

    /** */
    template <class Type>  CollectionNode<Type>& getCollection (const std::string& name)
    {
        /** See http://stackoverflow.com/questions/15572415/expected-primary-expression-before-in-g-but-not-in-microsoft-compiler */
        CollectionNode<Type>* result = Factory::template createCollection<Type> (this, name, 0);
        _collections.push_back (result);
        result->use ();
        return *result;
    }

    /** */
    void remove ()
    {
        /** We remove the collections. */
        for (size_t i=0; i<_collections.size(); i++)  {  _collections[i]->remove ();  }

        /** We remove the partitions. */
        for (size_t i=0; i<_partitions.size(); i++)   { _partitions[i]->remove(); }

        /** We remove the children groups. */
        for (size_t i=0; i<_groups.size(); i++)       { _groups[i]->remove(); }
    }

protected:

    std::vector<ICell*> _collections;
    std::vector<ICell*> _partitions;
    std::vector<Group<Factory>*> _groups;
};

/**********************************************************************
######      #     ######   #######  ###  #######  ###  #######  #     #
#     #    # #    #     #     #      #      #      #   #     #  ##    #
#     #   #   #   #     #     #      #      #      #   #     #  # #   #
######   #     #  ######      #      #      #      #   #     #  #  #  #
#        #######  #   #       #      #      #      #   #     #  #   # #
#        #     #  #    #      #      #      #      #   #     #  #    ##
#        #     #  #     #     #     ###     #     ###  #######  #     #
**********************************************************************/

/** \brief Define a set of Collection instances having the same type.
 *
 * The Partition class groups several Collections that share the same kind
 * of items.
 *
 * It is defined as a subclass of Group.
 */
template<typename Factory, typename Type>
class Partition : public Group<Factory>
{
public:

    /** Constructor.
     * \param[in] parent : the parent node
     * \param[in] id : the identifier of the instance to be created
     * \param[in] nbCollections : number of collections for this partition
     */
    Partition (dp::ICell* parent, const std::string& id, size_t nbCollections)
        : Group<Factory> (parent, id), _typedCollections(nbCollections), _synchro(0)
    {
        /** We create a synchronizer to be shared by the collections. */
        _synchro = system::impl::System::thread().newSynchronizer();

        /** We want to instantiate the wanted number of collections. */
        for (size_t i=0; i<_typedCollections.size(); i++)
        {
            /** We define the name of the current partition as a mere number. */
            std::stringstream ss;  ss << i;

            /** See http://stackoverflow.com/questions/15572415/expected-primary-expression-before-in-g-but-not-in-microsoft-compiler */
            CollectionNode<Type>* result = Factory::template createCollection<Type> (this, ss.str(), _synchro);

            /** We add the collection node to the dedicated vector and take a token for it. */
            (_typedCollections [i] = result)->use ();
        }
    }

    /** Destructor. */
    ~Partition ()
    {
        /** We release the token for each collection node. */
        for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->forget ();  }

        delete _synchro;
    }

    /** Return the number of collections for this partition.
     * \return the number of collections. */
    size_t size()  const  { return _typedCollections.size(); }

    /** Get the ith collection
     * \param[in] idx : index of the collection to be retrieved
     * \return the wanted collection.
     */
    Collection<Type>& operator[] (size_t idx)  {  return  * _typedCollections[idx]->getRef();  }

    /** Flush the whole partition (ie flush each collection). */
    void flush ()
    {
        for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->flush();  }
    }

    /** Remove physically the partition (ie. remove each collection). */
    void remove ()
    {
        /** We remove each collection of this partition. */
        for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->remove ();  }

        /** We call the remove of the parent class. */
        Group<Factory>::remove ();
    }

protected:

    std::vector <CollectionNode<Type>* > _typedCollections;
    system::ISynchronizer* _synchro;
};

/********************************************************************************/

/** \brief Cache (with potential synchronization) of a Partition.
 *
 * This class implements a cache for a Partition instance:
 *  -> an 'insert' is first cached in memory; when the cache is full, all the items are inserted in the
 *  real partition
 *  -> flushing the cache in the real partition is protected with a synchronizer
 *
 *  A typical use case is to create several PartitionCache referring the same Partition, and use them
 *  independently in different threads => in each thread, we will work on the local (cached) partition
 *  like a normal partition (ie. use 'insert'), but without taking care to the concurrent accesses to
 *  the referred partition (it is checked by the PartitionCache class).
 */
template<typename Factory, typename Type>
class PartitionCache
{
public:

    /** Constructor */
    PartitionCache (Partition<Factory,Type>& ref, size_t nbItemsCache, system::ISynchronizer* synchro=0)
        :  _ref(ref), _nbItemsCache(nbItemsCache), _synchro(synchro), _cachedCollections(ref.size())
    {
        /** We create the partition files. */
        for (size_t i=0; i<_cachedCollections.size(); i++)
        {
            _cachedCollections[i] = new CollectionCache<Type> (ref[i], nbItemsCache, synchro);
            _cachedCollections[i]->use ();
        }
    }

    /** Copy constructor. */
    PartitionCache (const PartitionCache<Factory, Type>& p)
        : _ref(p._ref), _nbItemsCache(p._nbItemsCache), _synchro(p._synchro), _cachedCollections(p.size())
    {
        /** We create the partition files. */
        for (size_t i=0; i<_cachedCollections.size(); i++)
        {
            PartitionCache<Factory,Type>& pp = (PartitionCache<Factory,Type>&)p;
            _cachedCollections[i] = new CollectionCache<Type> (pp[i].getRef(), p._nbItemsCache, p._synchro);
            _cachedCollections[i]->use ();
        }
    }

    /** Destructor. */
    ~PartitionCache ()
    {
        flush ();
        for (size_t i=0; i<_cachedCollections.size(); i++)  {  _cachedCollections[i]->forget ();  }
    }

    /** Return the number of collections for this partition.
     * \return the number of collections. */
    size_t size() const { return _cachedCollections.size(); }

    /** Get the ith collection
     * \param[in] idx : index of the collection to be retrieved
     * \return the wanted collection.
     */
    CollectionCache<Type>& operator[] (size_t idx)  {  return * _cachedCollections[idx];  }

    /** Flush the whole partition (ie flush each collection). */
    void flush ()   {  for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->flush();  }  }

    /** Remove physically the partition (ie. remove each collection). */
    void remove ()  {  for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->remove ();  } }

private:
    Partition<Factory,Type>& _ref;
    size_t                     _nbItemsCache;
    system::ISynchronizer*     _synchro;

    std::vector <CollectionCache<Type>* > _cachedCollections;
};

/********************************************************************************
        ######   ######   #######  ######   #     #   #####   #######
        #     #  #     #  #     #  #     #  #     #  #     #     #
        #     #  #     #  #     #  #     #  #     #  #           #
        ######   ######   #     #  #     #  #     #  #           #
        #        #   #    #     #  #     #  #     #  #           #
        #        #    #   #     #  #     #  #     #  #     #     #
        #        #     #  #######  ######    #####    #####      #
********************************************************************************/

/** \brief Product class
 *
 * The Product class is the entry point for managing collections and groups.
 *
 * It delegates all the actions to a root group (retrievable through an operator overload).
 *
 * Such a product is supposed to gather several information sets in a single environment,
 * with possible hierarchical composition.
 *
 * It is a template class: one should provide the actual type of the collections containers.
 * Possible template types could be designed for:
 *   - classical file system
 *   - HDF5 files
 *   - memory
 */
template<class Factory>
class Product : public dp::impl::Cell
{
public:

    /** Constructor.
     * \param[in] name : name of the product.
     * \param[in] autoRemove : tells whether the product has to be physically deleted when this object is deleted. */
    Product (const std::string& name, bool autoRemove=false)
        : dp::impl::Cell(0, ""),  _root(this, ""), _autoRemove(autoRemove)  {}

    /** Facility for retrieving the root group.
     * \return the root group. */
    Group<Factory>& operator() (const std::string name="")
    {
        if (name.empty())  { return _root; }
        else               { return _root.getGroup (name);  }
    }

    /** Remove physically the product. */
    virtual void remove ()  {  _root.remove(); }

protected:

    /** Root group. */
    Group<Factory> _root;

    bool _autoRemove;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_PRODUCT_HPP_ */
