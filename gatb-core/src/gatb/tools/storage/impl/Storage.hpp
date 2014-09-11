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

/** \file Storage.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Storage interface
 *
 *  This file holds interfaces related to the Collection, Group, Partition, Factory interfaces.
 *
 *  I believe this deals with both on-disk and in-memory structures.. (from discussions with Erwan)
 *  So if you're looking at this file to figure out whether a Collection/Group/Partition is stored on
 *  disk or on memory, look elsewhere!
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_PRODUCT_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_PRODUCT_HPP_

/********************************************************************************/

#include <gatb/tools/storage/impl/Cell.hpp>

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/CollectionAbstract.hpp>
#include <gatb/tools/collections/impl/CollectionCache.hpp>

#include <gatb/tools/misc/api/IProperty.hpp>

#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <map>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
/** File system storage for collections. */
namespace storage   {
namespace impl      {
/********************************************************************************/

/** */
enum StorageMode_e { STORAGE_FILE, STORAGE_HDF5, STORAGE_GZFILE, STORAGE_COMPRESSED_FILE };

/********************************************************************************/

/** Forward declarations. */
class Group;
template <typename Type>  class Partition;
template <typename Type>  class CollectionNode;

class StorageFactory;

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
class CollectionNode : public Cell, public collections::impl::CollectionAbstract<Item>
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
    CollectionNode (StorageFactory* factory, ICell* parent, const std::string& id, collections::Collection<Item>* ref);

    /** Destructor. */
    virtual ~CollectionNode();

    /** \copydoc dp::ICell::remove */
    void remove ();

    /** */
    void addProperty (const std::string& key, const std::string value);

    /** */
    std::string getProperty (const std::string& key);

    /**  */
    collections::Collection<Item>* getRef ();

private:

    StorageFactory* _factory;

    collections::Collection<Item>* _ref;
    void setRef (collections::Collection<Item>* ref)  { SP_SETATTR(ref); }
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
class Group : public Cell
{
public:

    /** Constructor. */
    Group (StorageFactory* factory, ICell* parent, const std::string& name);

    /** Destructor. */
    ~Group();

    /** */
    Group& getGroup (const std::string& name);

    /** */
    template <class Type>  Partition<Type>& getPartition (const std::string& name, size_t nb);

    /** */
    template <class Type>  CollectionNode<Type>& getCollection (const std::string& name);

    /** */
    void remove ();

    /** */
    virtual void addProperty (const std::string& key, const std::string value) { /*throw system::ExceptionNotImplemented ();*/ }

    /** */
    virtual std::string getProperty (const std::string& key)  { return "?"; /*throw system::ExceptionNotImplemented ();*/  }

protected:

    StorageFactory* _factory;
    std::vector<ICell*> _collections;
    std::vector<ICell*> _partitions;
    std::vector<Group*> _groups;
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
template<typename Type>
class Partition : public Group
{
public:

    /** Constructor.
     * \param[in] parent : the parent node
     * \param[in] id : the identifier of the instance to be created
     * \param[in] nbCollections : number of collections for this partition
     */
    Partition (StorageFactory* factory, ICell* parent, const std::string& id, size_t nbCollections);

    /** Destructor. */
    ~Partition ();

    /** Return the number of collections for this partition.
     * \return the number of collections. */
    size_t size()  const;

    /** Get the ith collection
     * \param[in] idx : index of the collection to be retrieved
     * \return the wanted collection.
     */
    collections::Collection<Type>& operator[] (size_t idx);

    /** Flush the whole partition (ie flush each collection). */
    void flush ();

    /** Remove physically the partition (ie. remove each collection). */
    void remove ();

protected:

    StorageFactory* _factory;
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
template<typename Type>
class PartitionCache
{
public:

    /** Constructor */
    PartitionCache (Partition<Type>& ref, size_t nbItemsCache, system::ISynchronizer* synchro=0);

    /** Copy constructor. */
    PartitionCache (const PartitionCache<Type>& p);

    /** Destructor. */
    ~PartitionCache ();

    /** Return the number of collections for this partition.
     * \return the number of collections. */
    size_t size() const;

    /** Get the ith collection
     * \param[in] idx : index of the collection to be retrieved
     * \return the wanted collection.
     */
    collections::impl::CollectionCache<Type>& operator[] (size_t idx);

    /** Flush the whole partition (ie flush each collection). */
    void flush ();

    /** Remove physically the partition (ie. remove each collection). */
    void remove ();

protected:
    Partition<Type>& _ref;
    size_t                     _nbItemsCache;
    system::ISynchronizer*     _synchro;
    std::vector <system::ISynchronizer*> _synchros;
    std::vector <collections::impl::CollectionCache<Type>* > _cachedCollections;
};

    
    
template<typename Type>
class PartitionCacheSorted
{
public:
    
    /** Constructor */
    PartitionCacheSorted (Partition<Type>& ref, size_t nbItemsCache, u_int32_t max_memory, system::ISynchronizer* synchro=0);
    
    /** Copy constructor. */
    PartitionCacheSorted (const PartitionCacheSorted<Type>& p);
    
    /** Destructor. */
    ~PartitionCacheSorted ();
    
    /** Return the number of collections for this partition.
     * \return the number of collections. */
    size_t size() const;
    
    /** Get the ith collection
     * \param[in] idx : index of the collection to be retrieved
     * \return the wanted collection.
     */
    collections::impl::CollectionCacheSorted<Type>& operator[] (size_t idx);
    
    /** Flush the whole partition (ie flush each collection). */
    void flush ();
    
    /** Remove physically the partition (ie. remove each collection). */
    void remove ();
    
protected:
    Partition<Type>& _ref;
    size_t                     _nbItemsCache;
    system::ISynchronizer*     _synchro;
    size_t _sharedBuffersSize;
    u_int32_t _max_memory;
    //only the model "parent" will have them, objects created with copy construc wont
    std::vector <system::ISynchronizer*> _synchros;
    std::vector <system::ISynchronizer*> _outsynchros;
    std::vector <Type*> _sharedBuffers;
    std::vector <size_t> _idxShared;
    int  _numthread;
    int  _nbliving;
    bool _own_synchros;
    bool _own_outsynchros;


    std::vector <collections::impl::CollectionCacheSorted<Type>* > _cachedCollections;
};
    

    
/********************************************************************************
        
           S t o r a g e

********************************************************************************/

/** \brief Storage class
 *
 * The Storage class is the entry point for managing collections and groups.
 *
 * It delegates all the actions to a root group (retrievable through an operator overload).
 *
 * Such a storage is supposed to gather several information sets in a single environment,
 * with possible hierarchical composition.
 *
 * It is a template class: one should provide the actual type of the collections containers.
 * Possible template types could be designed for:
 *   - classical file system
 *   - HDF5 files
 *   - memory
 */
class Storage : public Cell
{
public:

    /** Constructor.
     * \param[in] name : name of the storage.
     * \param[in] autoRemove : tells whether the storage has to be physically deleted when this object is deleted. */
    Storage (StorageMode_e mode, const std::string& name, bool autoRemove=false);

    /** Destructor */
    ~Storage ();

    Group& root () { return *getRoot(); }

    /** Facility for retrieving the root group.
     * \return the root group. */
    Group& getGroup (const std::string name) { return this->operator() (name); }

    /** Facility for retrieving the root group.
     * \return the root group. */
    Group& operator() (const std::string name="");

    /** Remove physically the storage. */
    virtual void remove ();

    /** */
    StorageFactory* getFactory() const { return _factory; }

protected:

    StorageFactory* _factory;
    void setFactory (StorageFactory* factory);

    /** Root group. */
    Group* _root;
    Group* getRoot ();
    void setRoot (Group* root)  { SP_SETATTR(root); }

    bool _autoRemove;
};

/********************************************************************************
                #######     #      #####   #######  #######  ######   #     #
                #          # #    #     #     #     #     #  #     #   #   #
                #         #   #   #           #     #     #  #     #    # #
   storage -    #####    #     #  #           #     #     #  ######      #
                #        #######  #           #     #     #  #   #       #
                #        #     #  #     #     #     #     #  #    #      #
                #        #     #   #####      #     #######  #     #     #
********************************************************************************/

class StorageFactory : public system::SmartPointer
{
public:

    /** Constructor */
    StorageFactory (StorageMode_e mode) : _mode(mode)  {}

    /** */
    Storage* create (const std::string& name, bool deleteIfExist, bool autoRemove);

    /** */
    bool exists (const std::string& name);

    /** */
    Storage* load (const std::string& name) { return create (name, false, false); }

    /** */
    Group* createGroup (ICell* parent, const std::string& name);

    /** */
    template<typename Type>
    Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb);

    /** */
    template<typename Type>
    CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro);

private:

    StorageMode_e _mode;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

/********************************************************************************/
/** WE INCLUDE THE 'IMPLEMENTATION' of the templates. */
#include <gatb/tools/storage/impl/Storage.tpp>
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_PRODUCT_HPP_ */
