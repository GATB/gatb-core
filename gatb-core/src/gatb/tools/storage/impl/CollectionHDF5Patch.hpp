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

/** \file CollectionHDF5Patch.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_HDF5_PATCH_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_HDF5_PATCH_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/CollectionAbstract.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>
#include <vector>
#include <stdarg.h>
#include <hdf5/hdf5.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

/* This class holds information to be shared by classes BagHDF5Patch and IterableHDF5Patch.
 *
 * IMPORTANT: HDF5 may suffer from kind of resources leaks when one uses H5Dread
 * several times; if no H5Dclose is done, some internal memory may be kept, looking
 * like memory leaks. See also:
 * http://stackoverflow.com/questions/18466691/excessive-memory-use-with-hdf5-h5dread
 *
 * Now, we use H5Dclose when the collection apparently doesn't look to be accessed
 * any more. If the datasetId is needed again, a request to get it again is needed.
 */
template <class Item> struct  CollectionDataHDF5Patch : public system::SmartPointer
{
//public:

    /** */
    CollectionDataHDF5Patch (hid_t fileId, const std::string& filename, system::ISynchronizer* synchro, int compress)
     : _fileId(fileId), _datasetId(0), _typeId(0), _nbItems(0), _name(filename), _synchro(synchro), _nbCalls(0), _compress(compress)
    {
        /** We get the HDF5 type of the item. */
        bool isCompound=false;
        _typeId = Item().hdf5(isCompound);

        /** We get the number of items in the data set.
         * NOTE: calling the (lazy) accessor 'getDatasetId' should also configure the dataset id. */
        hid_t filespaceId = H5Dget_space (this->getDatasetId());
        H5Sget_simple_extent_dims (filespaceId, &_nbItems, NULL);
        H5Sclose (filespaceId);
    }

    /** Destructor */
    ~CollectionDataHDF5Patch ()
    {
        herr_t status = H5Tclose (_typeId);
        if (status < 0)  { 
            std::cout << "HDF5 error (H5Tclose), status " <<  status << std::endl; exit(1); /* used to be an exception, but recent gcc's complain when in destructor*/  }

        if (_datasetId != 0) {  H5Dclose (_datasetId);   _datasetId=0; }
    }

    /**  */
    void addProperty (const std::string& key, const std::string value)
    {
        system::LocalSynchronizer localsynchro (_synchro);

        hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

        hsize_t dims = 1;
        hid_t space_id = H5Screate_simple (1, &dims, NULL);

        /** We create the attribute. */
        hid_t attrId = H5Acreate2 (getDatasetId(), key.c_str(), datatype,  space_id, H5P_DEFAULT, H5P_DEFAULT);

        /** We write the data. */
        const char* array[] = { value.c_str() };
        H5Awrite (attrId, datatype, &array);

        /** We close resources. */
        H5Aclose (attrId);
        H5Tclose (datatype);
        H5Sclose (space_id);
    }

    /** \copydoc tools::collections::Collection::getProperty */
    std::string getProperty (const std::string& key)
    {
        system::LocalSynchronizer localsynchro (_synchro);

        std::string result;

        herr_t status;

        hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

        hid_t attrId = H5Aopen (getDatasetId(), key.c_str(), H5P_DEFAULT);

        hid_t space_id = H5Aget_space (attrId);

        hsize_t dims = 1;
        H5Sget_simple_extent_dims (space_id, &dims, NULL);
        char** rdata = (char **) MALLOC (dims * sizeof (char *));

        status = H5Aread (attrId, datatype, rdata);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Aread), status %d", status);  }

        /** We set the result. */
        result.assign (rdata[0]);

        /** We release buffers. */
        status = H5Dvlen_reclaim (datatype, space_id, H5P_DEFAULT, rdata);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dvlen_reclaim), status %d", status);  }
        FREE (rdata);

        /** We close resources. */
        H5Aclose (attrId);
        H5Tclose (datatype);
        H5Sclose (space_id);

        return result;
    }

//protected:

    /** */
    hid_t retrieveDatasetId ()
    {
        hid_t result = 0;
        herr_t status;

        /** We look whether the object exists or not. */
        htri_t doesExist = H5Lexists (_fileId, _name.c_str(), H5P_DEFAULT);

        if (doesExist > 0)
        {
            result = H5Dopen2 (_fileId, _name.c_str(), H5P_DEFAULT);
        }
        else
        {
            /** We pack the type. */
            hid_t actualType = H5Tcopy (_typeId);

            /** We create the dataspace. */
            hsize_t dims      = 0;
            hsize_t maxdims   = H5S_UNLIMITED;
            hid_t dataspaceId = H5Screate_simple (1, &dims, &maxdims);

            /* Modify dataset creation properties, i.e. enable chunking  */
            hsize_t chunk_dims = GATB_HDF5_NB_ITEMS_PER_BLOCK;
            hid_t propId = H5Pcreate     (H5P_DATASET_CREATE);

            if (_compress > 0)
            {
            	status = H5Pset_shuffle (propId);
            	status = H5Pset_deflate (propId, _compress);
            }

            status       = H5Pset_layout (propId, H5D_CHUNKED);
            status       = H5Pset_chunk  (propId, 1, &chunk_dims);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Pset_chunk), status %d", status);  }

            /** We create the dataset. */
            result = H5Dcreate2 (_fileId, _name.c_str(),  actualType, dataspaceId, H5P_DEFAULT, propId, H5P_DEFAULT);

            /** Some cleanup. */
            H5Pclose (propId);
            H5Sclose (dataspaceId);
            H5Tclose (actualType);
        }

        return result;
    }

    hid_t getDatasetId ()
    {
        if (this->_datasetId == 0)  {  this->_datasetId = this->retrieveDatasetId ();  }
        return this->_datasetId;
    }

    /** */
    void clean ()
    {
        if (_datasetId == 0)  { return; }
        {
            system::LocalSynchronizer localsynchro (_synchro);
            if (_datasetId != 0) {  H5Dclose (_datasetId);   _datasetId=0; }
        }
    }

    hid_t                   _fileId;
    hid_t                   _datasetId;
    hid_t                   _typeId;
    hsize_t                 _nbItems;
    std::string             _name;
    system::ISynchronizer*  _synchro;
    u_int64_t               _nbCalls;
    int                     _compress;

    void checkCleanup ()
    {
        static u_int64_t MASK = ( (1<<GATB_HDF5_CLEANUP_WORKAROUND) - 1);

        u_int64_t newValue = __sync_add_and_fetch (&_nbCalls, 1);

        /** We periodically clean up some HDF5 resources. */
        if ( (newValue & MASK) == MASK) {  clean();  }
    }
};

/********************************************************************************/

template <class Item> class BagHDF5Patch : public collections::Bag<Item>, public system::SmartPointer
{
public:

    /** Constructor */
    BagHDF5Patch (CollectionDataHDF5Patch<Item>* common)  : _common (0), _nbInserted(0)  { setCommon(common); }

    /** Destructor. */
    ~BagHDF5Patch ()  { setCommon(0); }

    /** Insert an item into the bag.
     * \param[in] item : the item to be inserted. */
    void insert (const Item& item) {  insert (&item, 1);  }

    void insert (const std::vector<Item>& items, size_t length=0)  {  insert (items.data(), length==0 ? items.size() : length); }

    /** Insert items into the bag.
     * \param[in] items : items to be inserted. */
    void insert (const Item* items, size_t length)
    {
        if (items==0 || length==0)  { return; }
        
        herr_t status = 0;

        {
            system::LocalSynchronizer localsynchro (_common->_synchro);

            /** We get the dataset id. */
            hid_t datasetId = _common->getDatasetId();

            /** Resize the memory dataspace to indicate the new size of our buffer. */
            hsize_t memDim = length;
            hid_t memspaceId = H5Screate_simple (1, &memDim, NULL);

            /** Extend dataset. */
            hsize_t newDim = _nbInserted + length;
            status = H5Dset_extent (datasetId, &newDim);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dset_extent), status %d", status);  }

            /** Select hyperslab on file dataset. */
            hid_t filespaceId = H5Dget_space(datasetId);
            hsize_t start = _nbInserted;
            hsize_t count = length;
            status = H5Sselect_hyperslab (filespaceId, H5S_SELECT_SET, &start, NULL, &count, NULL);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sselect_hyperslab), status %d", status);  }

            /** Append buffer to dataset */
            status = H5Dwrite (datasetId, _common->_typeId, memspaceId, filespaceId, H5P_DEFAULT, items);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dwrite), status %d", status);  }

            /** We increase the number of inserted items. */
            _nbInserted += length;

            __sync_fetch_and_add (&_common->_nbItems, length);

            /** Close resources. */
            status = H5Sclose (filespaceId);
            status = H5Sclose (memspaceId);
            if (status != 0)  { std::cout << "err H5Sclose" << std::endl; }
        }

        /** We periodically clean up some HDF5 resources. */
        _common->checkCleanup ();
    }

    /** */
    void flush ()  {  }

    CollectionDataHDF5Patch<Item>* _common;
    void setCommon (CollectionDataHDF5Patch<Item>* common)  { SP_SETATTR(common); }

private:

    u_int64_t _nbInserted;
};

/********************************************************************************/

template<typename Item> class HDF5IteratorPatch;

template <class Item> class IterableHDF5Patch : public collections::Iterable<Item>,  public system::SmartPointer
{
public:

    /** */
    IterableHDF5Patch (CollectionDataHDF5Patch<Item>* common)  : _common (0) { setCommon(common); }

    /** */
    ~IterableHDF5Patch ()  { setCommon(0);}

    /** */
    dp::Iterator<Item>* iterator ()  {  return new HDF5IteratorPatch<Item> (this);  }

    /** */
    int64_t getNbItems ()  {  return _common->_nbItems;  }

    /** */
    int64_t estimateNbItems ()  {  return getNbItems();  }

    /** */
    Item* getItems (Item*& buffer)
    {
        retrieveCache (buffer, 0, getNbItems());
        return buffer;
    }

    /** */
    size_t getItems (Item*& buffer, size_t start, size_t count)
    {
        return retrieveCache (buffer, start, count);
    }

private:

    CollectionDataHDF5Patch<Item>* _common;
    void setCommon (CollectionDataHDF5Patch<Item>* common)  { SP_SETATTR(common); }

    /** */
    u_int64_t retrieveCache (Item* data, hsize_t start, hsize_t count)
    {
        herr_t status = 0;

        if (start       > _common->_nbItems)  { return 0;                          }
        if (start+count > _common->_nbItems)  { count = _common->_nbItems - start; }

        /** We use a synchronizer instruction block.
         * NOTE !!!  the 'clean' method called after this block is also synchronized,
         * and therefore must not be in the same instruction block. */
        {
            system::LocalSynchronizer localsynchro (_common->_synchro);

            hid_t memspaceId = H5Screate_simple (1, &count, NULL);

            /** Select hyperslab on file dataset. */
            hid_t filespaceId = H5Dget_space(_common->getDatasetId());
            status = H5Sselect_hyperslab (filespaceId, H5S_SELECT_SET, &start, NULL, &count, NULL);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sselect_hyperslab), status %d", status);  }

            /** Read buffer from dataset */
            status = H5Dread (_common->getDatasetId(), _common->_typeId, memspaceId, filespaceId, H5P_DEFAULT, data);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dread), status %d", status);  }

            /** Close resources. */
            status = H5Sclose (filespaceId);
            status = H5Sclose (memspaceId);
            if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sclose), status %d", status);  }
        }

        /** We periodically clean up some HDF5 resources. */
        _common->checkCleanup ();

        return count;
    }

    template<typename U>
    friend class HDF5IteratorPatch;
};

/********************************************************************************/

template<typename Item>
class HDF5IteratorPatch : public dp::Iterator<Item>
{
public:

    /** */
    HDF5IteratorPatch ()  : _ref(0), _blockSize(0),
         _data(0), _dataSize(0), _dataIdx(0), _isDone (true),
         _nbRead(0), _total(0), _memspaceId(0)
    {}

    /** */
    HDF5IteratorPatch (const HDF5IteratorPatch& it)
    : _ref(it._ref), _blockSize(it._blockSize),
        _data(0), _dataSize(0), _dataIdx(0), _isDone (true),
        _nbRead(0), _total(0), _memspaceId(0)
    {
        _data = (Item*) MALLOC (_blockSize*sizeof(Item));
        memset (_data, 0, _blockSize*sizeof(Item));
        _total = _ref->_nbItems;
    }

    /** */
    HDF5IteratorPatch (IterableHDF5Patch<Item>* ref, size_t blockSize=GATB_HDF5_NB_ITEMS_PER_BLOCK)
        : _ref(ref), _blockSize(blockSize),
          _data(0), _dataSize(0), _dataIdx(0), _isDone (true),
          _nbRead(0), _total(0), _memspaceId(0)
    {
        _data = (Item*) MALLOC (_blockSize*sizeof(Item));
        memset (_data, 0, _blockSize*sizeof(Item));

        _total = _ref->_common->_nbItems;
    }

    HDF5IteratorPatch& operator= (const HDF5IteratorPatch& it)
    {
        if (this != &it)
        {
            _ref        = it._ref;
            _blockSize  = it._blockSize;
            _dataSize   = it._dataSize;
            _dataIdx    = it._dataIdx;
            _isDone     = it._isDone;
            _nbRead     = it._nbRead;
            _memspaceId = it._memspaceId;
            _total      = it._total;

            if (_data)  { FREE (_data); }
            _data = (Item*) MALLOC (_blockSize*sizeof(Item));
            memcpy (_data, it._data, _blockSize*sizeof(Item));
        }
        return *this;
    }

    /** */
    ~HDF5IteratorPatch()
    {
        if (_data)  { FREE (_data); }

        /** We clean the reference instance. */
        _ref->_common->clean();
    }

    void first()
    {
        _nbRead   = 0;
        _dataIdx  = 0;
        _dataSize = retrieveNextCache();
        _isDone   = _dataIdx >= _dataSize;

        if (!_isDone)  {  *this->_item = _data[_dataIdx];  }
    }

    void next()
    {
        _dataIdx++;
        _isDone = _dataIdx >= _dataSize;
        if (!_isDone)  {  *this->_item = _data[_dataIdx];  }
        else
        {
            _dataIdx  = 0;
            _dataSize = retrieveNextCache();
            _isDone = _dataIdx >= _dataSize;
            if (!_isDone)  {  *this->_item = _data[_dataIdx];  }
        }
    }

    bool isDone()
    {
        return _isDone;
    }

    Item& item ()  { return *this->_item; }

private:
    IterableHDF5Patch<Item>* _ref;
    size_t        _blockSize;

    Item*     _data;
    u_int64_t _dataSize;
    u_int64_t _dataIdx;

    bool _isDone;

    u_int64_t     _nbRead;
    u_int64_t     _total;

    hid_t _memspaceId;

    u_int64_t retrieveNextCache ()
    {
        if (_total <= _nbRead)  {  return 0;   }
        hsize_t nbToRead = std::min ((u_int64_t)_blockSize, _total - _nbRead);

        _nbRead += _ref->retrieveCache (_data, _nbRead, nbToRead);

        return nbToRead;
    }
};

/********************************************************************************/

template <class Item> class CollectionHDF5Patch : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionHDF5Patch (hid_t fileId, const std::string& name, system::ISynchronizer* synchro, int compressLevel)
        : collections::impl::CollectionAbstract<Item> (0,0), _common(0)
    {
        system::LocalSynchronizer localsynchro (synchro);

        CollectionDataHDF5Patch<Item>* common = new CollectionDataHDF5Patch<Item> (fileId, name, synchro, compressLevel);

        /** We create the bag and the iterable instances. */
        this->setBag      (new BagHDF5Patch<Item>      (common));
        this->setIterable (new IterableHDF5Patch<Item> (common));
    }

    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {}

    /** \copydoc tools::collections::Collection::addProperty */
    void addProperty (const std::string& key, const std::string value)
    {
        BagHDF5Patch<Item>* theBag = dynamic_cast<BagHDF5Patch<Item>*> (this->bag());
        if (theBag != 0)  {  theBag->_common->addProperty (key, value);  }
    }

    /** \copydoc tools::collections::Collection::getProperty */
    std::string getProperty (const std::string& key)
    {
        std::string result;

        BagHDF5Patch<Item>* theBag = dynamic_cast<BagHDF5Patch<Item>*> (this->bag());
        if (theBag != 0)  {  result = theBag->_common->getProperty (key);  }

        return result;
    }

private:
    CollectionDataHDF5Patch<Item>* _common;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_HDF5_PATCH_HPP_ */
