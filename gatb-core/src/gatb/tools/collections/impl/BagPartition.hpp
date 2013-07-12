/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file BagPartition.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bag partition implementation for file
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_PARTITION_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_PARTITION_HPP_

/********************************************************************************/

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>
#include <vector>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

template<typename Item> class BagFilePartition
{
public:
    BagFilePartition (size_t nbPartitions, const std::string& format) : _partitions(nbPartitions), _uriFormat(format)
    {
        /** We delete physically the partition files. */
        for (size_t i=0; i<_partitions.size(); i++)  {  system::impl::System::file().remove (getFilename(i));  }

        /** We create the partition files. */
        for (size_t i=0; i<_partitions.size(); i++)
        {
            _partitions[i] = new BagFile<Item> (getFilename(i));
            _partitions[i]->use ();
        }
    }

    /** */
    virtual ~BagFilePartition ()
    {
        /** We logically delete the partition files. */
        for (size_t i=0; i<_partitions.size(); i++)
        {
            _partitions[i]->forget();
        }
    }

    /** */
    Bag<Item>*& operator[] (size_t idx)  { return _partitions[idx]; }

    /** */
    size_t size () const  { return _partitions.size(); }

private:

    std::string getFilename (size_t idx)
    {
        char filename[128];  snprintf (filename, sizeof(filename), _uriFormat.c_str(), idx);
        return std::string(filename);
    }

    std::vector<Bag<Item>*> _partitions;
    std::string             _uriFormat;
};

/********************************************************************************/

template<typename Item> class BagCachePartition
{
public:
    BagCachePartition (BagFilePartition<Item>& partition, system::ISynchronizer* synchro)
        : _ref(partition), cache(partition.size()), _synchro(synchro), _cacheNbItems(1<<12)
    {
        /** We create the partition files. */
        for (size_t i=0; i<cache.size(); i++)
        {
            cache[i] = createBag ((partition[i]), _cacheNbItems, synchro);
            cache[i]->use ();
        }
    }

    BagCachePartition (const BagCachePartition<Item>& p)
        : _ref(p._ref), cache(p.cache.size()), _synchro(p._synchro), _cacheNbItems(p._cacheNbItems)
    {
        /** We create the partition files. */
        for (size_t i=0; i<cache.size(); i++)
        {
            cache[i] = createBag ((p._ref[i]), _cacheNbItems, _synchro);
            cache[i]->use ();
        }
    }

    virtual ~BagCachePartition ()
    {
        for (size_t i=0; i<cache.size(); i++)
        {
            cache[i]->flush();
            cache[i]->forget ();
        }
    }

    /** */
    Bag<Item>*& operator[] (size_t idx)  { return cache[idx]; }

    /** */
    size_t size () const  { return cache.size(); }

private:

    virtual BagCache<Item>* createBag (Bag<Item>*& b, size_t cacheSize, system::ISynchronizer* synchro)
    {
        return new BagCache<Item> (b, cacheSize, synchro);
        //return new BagCacheSorted<Item> (b, cacheSize, synchro);  // don't use when using OAHash for DSK solid kmers filling.
    }

    BagFilePartition<Item>&  _ref;
    std::vector <Bag<Item>*>  cache;
    system::ISynchronizer*   _synchro;
    size_t                   _cacheNbItems;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_PARTITION_HPP_ */
