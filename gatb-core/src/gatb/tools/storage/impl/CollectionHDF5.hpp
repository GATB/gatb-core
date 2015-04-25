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

/** \file CollectionHDF5.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_HDF5_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_HDF5_HPP_

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

/** \brief Implementation of the Bag interface with a HDF5 file.
 *
 * This implementation writes Item objects in a HDF5 file.
 */
template <class Item> class BagHDF5 : public collections::Bag<Item>, public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] datasetId : HDF5 identifier of the dataset acting as a Bag.
     * \param[in] typeId : HDF5 type identifier for the Item type
     * \param[in] nbItems : number of items
     * \param[in] synchro : used to serialize concurrent read/write HDF5 operations.
     */
    BagHDF5 (hid_t datasetId, hid_t typeId, u_int64_t& nbItems, system::ISynchronizer* synchro)
        : _datasetId(datasetId), _typeId(typeId), _nbInserted(0), _nbItems(nbItems), _synchro(synchro)
    {
    }

    /** \copydoc collections::Bag::insert */
    void insert (const Item& item) {  insert (&item, 1);  }

    /** \copydoc collections::Bag::insert(const std::vector<Item>& items, size_t length) */
    void insert (const std::vector<Item>& items, size_t length=0)  {  insert (items.data(), length==0 ? items.size() : length); }

    /** \copydoc collections::Bag::insert(const Item* items, size_t length) */
    void insert (const Item* items, size_t length)
    {
        if (items==0 || length==0)  { return; }
        
        herr_t status = 0;

        system::LocalSynchronizer localsynchro (_synchro);

        /** Resize the memory dataspace to indicate the new size of our buffer. */
        hsize_t memDim = length;
        hid_t memspaceId = H5Screate_simple (1, &memDim, NULL);

        /** Extend dataset. */
        hsize_t newDim = _nbInserted + length;
        status = H5Dset_extent (_datasetId, &newDim);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dset_extent), status %d", status);  }

        /** Select hyperslab on file dataset. */
        hid_t filespaceId = H5Dget_space(_datasetId);
        hsize_t start = _nbInserted;
        hsize_t count = length;
        status = H5Sselect_hyperslab (filespaceId, H5S_SELECT_SET, &start, NULL, &count, NULL);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sselect_hyperslab), status %d", status);  }

        /** Append buffer to dataset */
        status = H5Dwrite (_datasetId, _typeId, memspaceId, filespaceId, H5P_DEFAULT, items);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dwrite), status %d", status);  }

        /** We increase the number of inserted items. */
        _nbInserted += length;

        __sync_fetch_and_add (&_nbItems, length);

        /** Close resources. */
        status = H5Sclose (filespaceId);
        status = H5Sclose (memspaceId);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sclose), status %d", status);  }

    }

    /** \copydoc collections::Bag::flush */
    void flush ()  {}

private:

    hid_t     _datasetId;
    hid_t     _typeId;
    u_int64_t& _nbItems;
    u_int64_t _nbInserted;
    system::ISynchronizer* _synchro;

};

/********************************************************************************/

template<typename Item> class HDF5Iterator;

template <class Item> class IterableHDF5 : public collections::Iterable<Item>, public system::SmartPointer
{
public:

    /** */
    IterableHDF5 (hid_t datasetId, hid_t typeId, u_int64_t& nbItems,     system::ISynchronizer* synchro)
        :  _datasetId(datasetId), _typeId(typeId), _nbItems(nbItems), _synchro(synchro)  {}

    /** */
    ~IterableHDF5 ()  {}

    /** */
    dp::Iterator<Item>* iterator ()  {  return new HDF5Iterator<Item> (this);  }

    /** */
    int64_t getNbItems ()
    {
        return _nbItems;
    }

    /** */
    int64_t estimateNbItems ()
    {
        return getNbItems();
    }

    /** */
    Item* getItems (Item*& buffer)
    {
        retrieveCache (buffer, 0, getNbItems());
        return buffer;
    }

    size_t getItems (Item*& buffer, size_t start, size_t count)  { return retrieveCache (buffer, start, count); }


private:

    hid_t _datasetId;
    hid_t _typeId;
    u_int64_t& _nbItems;
    system::ISynchronizer* _synchro;


    /** */
    u_int64_t retrieveCache (Item* data, hsize_t start, hsize_t count)
    {
        herr_t status = 0;

        if (start       > _nbItems)  { return 0;               }
        if (start+count > _nbItems)  { count = _nbItems-start; }

        system::LocalSynchronizer localsynchro (_synchro);

        hid_t memspaceId = H5Screate_simple (1, &count, NULL);

        /** Select hyperslab on file dataset. */
        hid_t filespaceId = H5Dget_space(_datasetId);
        status = H5Sselect_hyperslab (filespaceId, H5S_SELECT_SET, &start, NULL, &count, NULL);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sselect_hyperslab), status %d", status);  }

        /** Read buffer from dataset */
        status = H5Dread (_datasetId, _typeId, memspaceId, filespaceId, H5P_DEFAULT, data);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Dread), status %d", status);  }

        /** Close resources. */
        status = H5Sclose (filespaceId);
        status = H5Sclose (memspaceId);
        if (status < 0)  { throw gatb::core::system::Exception ("HDF5 error (H5Sclose), status %d", status);  }

        return count;
    }

    template<typename U>
    friend class HDF5Iterator;
};

/********************************************************************************/

/** */
template<typename Item>
class HDF5Iterator : public dp::Iterator<Item>
{
public:

    /** */
    HDF5Iterator ()  : _ref(0), _blockSize(0),
         _data(0), _dataSize(0), _dataIdx(0), _isDone (true),
         _nbRead(0), _memspaceId(0), _total(0)
    {}

    /** */
    HDF5Iterator (const HDF5Iterator& it)
    : _ref(it._ref), _blockSize(it._blockSize),
        _data(0), _dataSize(0), _dataIdx(0), _isDone (true),
        _nbRead(0), _memspaceId(0), _total(0)
    {
        _data = (Item*) MALLOC (_blockSize*sizeof(Item));
        memset (_data, 0, _blockSize*sizeof(Item));
        _total = _ref->_nbItems;
    }

    /** */
    HDF5Iterator (IterableHDF5<Item>* ref, size_t blockSize=4096)
        : _ref(ref), _blockSize(blockSize),
          _data(0), _dataSize(0), _dataIdx(0), _isDone (true),
          _nbRead(0), _memspaceId(0), _total(0)
    {
        _data = (Item*) MALLOC (_blockSize*sizeof(Item));
        memset (_data, 0, _blockSize*sizeof(Item));

        _total = _ref->_nbItems;
    }

    HDF5Iterator& operator= (const HDF5Iterator& it)
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
    ~HDF5Iterator()
    {
        if (_data)  { FREE (_data); }
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
    IterableHDF5<Item>* _ref;
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

/** \brief Collection interface
 */
template <class Item> class CollectionHDF5 : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionHDF5 (hid_t fileId, const std::string& filename, system::ISynchronizer* synchro)
        : collections::impl::CollectionAbstract<Item> (0,0), _datasetId(0), _typeId(0), _nbItems(0)
    {
        herr_t status;

        system::LocalSynchronizer localsynchro (synchro);

        /** We get the HDF5 type of the item. */
        bool isCompound=false;
        _typeId = Item().hdf5(isCompound);

        /** We pack the type. */
        hid_t actualType = H5Tcopy (_typeId);
        //if (isCompound)  {  status = H5Tpack(actualType);  }

        std::string actualName = filename;

        /** We look whether the object exists or not. */
        htri_t doesExist = H5Lexists (fileId, actualName.c_str(), H5P_DEFAULT);

        if (doesExist > 0)
        {
            _datasetId = H5Dopen2 (fileId, actualName.c_str(), H5P_DEFAULT);

            hid_t filespaceId = H5Dget_space (_datasetId);

            hsize_t dims;
            H5Sget_simple_extent_dims (filespaceId, &dims, NULL);

            _nbItems = dims;
            status = H5Sclose (filespaceId);
        }
        else
        {
            /** We create the dataspace. */
            hsize_t dims      = 0;
            hsize_t maxdims   = H5S_UNLIMITED;
            hid_t dataspaceId = H5Screate_simple (1, &dims, &maxdims);

            /* Modify dataset creation properties, i.e. enable chunking  */
            hsize_t chunk_dims = 4096;
            hid_t propId = H5Pcreate     (H5P_DATASET_CREATE);
            status       = H5Pset_layout (propId, H5D_CHUNKED);
            status       = H5Pset_chunk  (propId, 1, &chunk_dims);

            /** We create the dataset. */
            _datasetId = H5Dcreate2 (fileId, filename.c_str(),  actualType, dataspaceId, H5P_DEFAULT, propId, H5P_DEFAULT);

            /** Some cleanup. */
            H5Pclose (propId);
            H5Sclose (dataspaceId);
            H5Tclose (actualType);
        }

        /** We create the bag and the iterable instances. */
        this->setBag      (new BagHDF5<Item>      (_datasetId, _typeId, _nbItems, synchro));
        this->setIterable (new IterableHDF5<Item> (_datasetId, _typeId, _nbItems, synchro));

    }

    /** Destructor. */
    virtual ~CollectionHDF5()
    {
        herr_t status;
        status = H5Dclose (_datasetId);
        status = H5Tclose (_typeId);
    }

    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {}

    /** \copydoc tools::collections::Collection::addProperty */
    void addProperty (const std::string& key, const std::string value)
    {
        hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

        hsize_t dims = 1;
        hid_t space_id = H5Screate_simple (1, &dims, NULL);

        /** We create the attribute. */
        hid_t attrId = H5Acreate2 (_datasetId, key.c_str(), datatype,  space_id, H5P_DEFAULT, H5P_DEFAULT);

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
        std::string result;
        herr_t status;

        hid_t datatype = H5Tcopy (H5T_C_S1);  H5Tset_size (datatype, H5T_VARIABLE);

        hid_t attrId = H5Aopen (_datasetId, key.c_str(), H5P_DEFAULT);

        hid_t space_id = H5Aget_space (attrId);

        hsize_t dims = 1;
        H5Sget_simple_extent_dims (space_id, &dims, NULL);
        char** rdata = (char **) MALLOC (dims * sizeof (char *));

        status = H5Aread (attrId, datatype, rdata);

        /** We set the result. */
        result.assign (rdata[0]);

        /** We release buffers. */
        status = H5Dvlen_reclaim (datatype, space_id, H5P_DEFAULT, rdata);
        free (rdata);

        /** We close resources. */
        H5Aclose (attrId);
        H5Tclose (datatype);
        H5Sclose (space_id);

        return result;
    }

private:

    hid_t _datasetId;
    hid_t _typeId;
    u_int64_t _nbItems;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_HDF5_HPP_ */
