/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file BagCache.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bag implementation as a cache for a referred Bag
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/system/impl/System.hpp>
#include <vector>
#include <algorithm>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/** Implementation of API from collections */
namespace impl          {
/********************************************************************************/

/** \brief Bag implementation as a cache to a referred Bag instance
 */
template <typename Item> class BagCache : public Bag<Item>
{
public:

    /** Constructor. */
    BagCache (Bag<Item>& ref, size_t cacheSize, system::ISynchronizer* synchro=0)
        : _ref(ref), _nbMax(cacheSize), _synchro(synchro), _items (cacheSize), _idx(0)  {}

    /** Destructor. */
    ~BagCache ()  {}

    /**  \copydoc Bag::insert */
    void insert (const Item& item)
    {
        if (_idx+1 > _nbMax)
        {
            if (_synchro)  {  _synchro->lock();    }
            flushCache ();
            if (_synchro)  {  _synchro->unlock();  }
        }

        _items[_idx++] = item;
    }

    /**  \copydoc Bag::flush */
    void flush ()
    {
        if (_synchro)  {  _synchro->lock();    }
        flushCache ();
        _ref.flush();
        if (_synchro)  {  _synchro->unlock();  }
    }

private:
    Bag<Item>&             _ref;
    system::ISynchronizer* _synchro;
    std::vector<Item>      _items;
    size_t                 _nbMax;
    size_t                 _idx;

    void flushCache ()
    {
        if (_idx > 0)  { _ref.insert (_items, _idx);  }
        _idx = 0;
    }
};

/********************************************************************************/

/** \brief Bag implementation as a cache to a referred Bag instance.
 *
 * The cache is sorted before sent to the reference
 */
template <typename Item> class BagCacheSorted : public Bag<Item>
{
public:

    /** Constructor. */
    BagCacheSorted (Bag<Item>& ref, size_t cacheSize, system::ISynchronizer* synchro=0)
: _ref(ref), _nbMax(cacheSize), _synchro(synchro), _items (cacheSize), _idx(0)  {}

    /** Destructor. */
    ~BagCacheSorted ()  {}

    /**  \copydoc Bag::insert */
    void insert (const Item& item)
    {
        if (_idx+1 > _nbMax)
        {
            std::sort (_items.begin (), _items.end ());

            if (_synchro)  {  _synchro->lock();      }
            flushCache ();
            if (_synchro)  {  _synchro->unlock();    }
        }

        _items[_idx++] = item;
    }

    /**  \copydoc Bag::flush */
    void flush ()
    {
        _items.resize (_idx);
        std::sort (_items.begin (), _items.end ());

        if (_synchro)  {  _synchro->lock();    }
        flushCache ();
        _ref.flush();
        if (_synchro)  {  _synchro->unlock();  }
    }

private:
    Bag<Item>&             _ref;
    system::ISynchronizer* _synchro;
    std::vector<Item>      _items;
    size_t                 _nbMax;
    size_t                 _idx;

    void flushCache ()
    {
        _ref.insert (_items, _idx);
        _idx = 0;
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_ */
