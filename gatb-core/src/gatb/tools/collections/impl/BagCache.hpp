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
    BagCache (Bag<Item>* ref, size_t cacheSize, system::ISynchronizer* synchro=0)
        : _ref(0), _nbMax(cacheSize), _synchro(synchro), _items (cacheSize), _idx(0)
    {
        setRef(ref);
    }

    BagCache (const BagCache<Item>& b)
        : _ref(0), _nbMax(b._nbMax), _synchro(b._synchro), _items (b._nbMax), _idx(0)

    {
        setRef (b._ref);
    }


    /** Destructor. */
    ~BagCache ()
    {
        setRef(0);
    }

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
        _ref->flush();
        if (_synchro)  {  _synchro->unlock();  }
    }

protected:

    Bag<Item>*             _ref;
    void setRef (Bag<Item>* ref)  { SP_SETATTR(ref); }

    system::ISynchronizer* _synchro;
    std::vector<Item>      _items;
    size_t                 _nbMax;
    size_t                 _idx;

    void flushCache ()
    {
        if (_idx > 0)  { _ref->insert (_items, _idx);  }
        _idx = 0;
    }
};

/********************************************************************************/

/** \brief Bag implementation as a cache to a referred Bag instance.
 *
 * The cache is sorted before sent to the reference
 */
template <typename Item> class BagCacheSorted : public BagCache<Item>
{
public:

    /** Constructor. */
    BagCacheSorted (Bag<Item>* ref, size_t cacheSize, system::ISynchronizer* synchro=0)
        : BagCache<Item>(ref, cacheSize, synchro)  { }

    /**  \copydoc Bag::insert */
    void insert (const Item& item)
    {
        if (this->_idx+1 > this->_nbMax)
        {
            std::sort (this->_items.begin (), this->_items.end ());

            if (this->_synchro)  {  this->_synchro->lock();      }
            this->flushCache ();
            if (this->_synchro)  {  this->_synchro->unlock();    }
        }

        this->_items[this->_idx++] = item;
    }

    /**  \copydoc Bag::flush */
    void flush ()
    {
        this->_items.resize (this->_idx);
        std::sort (this->_items.begin (), this->_items.end ());

        if (this->_synchro)  {  this->_synchro->lock();    }
        this->flushCache ();
        this->_ref->flush();
        if (this->_synchro)  {  this->_synchro->unlock();  }
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_ */
