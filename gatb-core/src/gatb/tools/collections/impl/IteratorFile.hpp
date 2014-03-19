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

/** \file IteratorFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterator implementation for file
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>
#include <vector>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/** Implementation of API from collections */
namespace impl          {
/********************************************************************************/

#define BUFFER_SIZE (128*1024)

/** \brief Iterator implementation for file
 */
template <class Item> class IteratorFile : public dp::Iterator<Item>
{
public:

    /** Constructor. */
    IteratorFile () : _file(0), _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(0), _isDone(true) {}

    IteratorFile (const IteratorFile& it):
        _filename(it._filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(it._cacheItemsNb), _isDone(true)
    {
        _file    = system::impl::System::file().newFile (_filename, "rb");
        _buffer  = (Item*) malloc (sizeof(Item) * _cacheItemsNb);
    }

    /** Constructor. */
    IteratorFile (const std::string& filename, size_t cacheItemsNb=10000) :
        _filename(filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(cacheItemsNb), _isDone(true)

    {
        _file    = system::impl::System::file().newFile (filename, "rb");
        _buffer  = (Item*) malloc (sizeof(Item) * _cacheItemsNb);
    }

    /** Destructor. */
    ~IteratorFile ()
    {
        if (_file)  { delete _file;  }
        if (_buffer) { free (_buffer); }
    }

    /** Affectation. */
    IteratorFile& operator= (const IteratorFile& it)
    {
        if (this != &it)
        {
            if (_file)    { delete _file; }
            if (_buffer)  { free(_buffer); }

            _filename     = it._filename;
            _cpt_buffer   = it._cpt_buffer;
            _idx          = it._idx;
            _cacheItemsNb = it._cacheItemsNb;
            _isDone       = it._isDone;

            _file    = system::impl::System::file().newFile (it._filename, "rb");
            _buffer  = (Item*) malloc (sizeof(Item) * it._cacheItemsNb);
        }
        return *this;
    }

    /** \copydoc dp::Iterator::first */
    void first()
    {
        _file->seeko (0, SEEK_SET);
        _cpt_buffer = 0;
        _idx        = 0;
        _isDone     = false;
        next ();
    }

    /** \copydoc dp::Iterator::next */
    void next()
    {
        if (_cpt_buffer==0)
        {
            _idx = 0;
            _cpt_buffer = _file->fread (_buffer, sizeof(Item), _cacheItemsNb);
            if (_cpt_buffer==0)  { _isDone = true;  return; }
        }

        *(this->_item) =  _buffer[_idx];
        _cpt_buffer --;
        _idx ++;
    }

    /** \copydoc dp::Iterator::isDone */
    bool isDone()  { return _isDone; }

    /** \copydoc dp::Iterator::item */
    Item& item ()  { return *(this->_item); }

    /** */
    size_t fill (std::vector<Item>& vec, size_t len=0)
    {
        if (len==0)  { len = vec.size(); }
        return _file->fread (vec.data(), sizeof(Item), len);
    }

private:
    std::string     _filename;
    system::IFile*  _file;
    Item*           _buffer;
    int             _cpt_buffer;
    int             _idx;
    size_t          _cacheItemsNb;
    bool            _isDone;
};

/********************************************************************************/

template <class Item> class IterableFile : public tools::collections::Iterable<Item>, public virtual system::SmartPointer
{
public:

    /** */
    IterableFile (const std::string& filename, size_t cacheItemsNb=10000)
        :   _filename(filename), _cacheItemsNb (cacheItemsNb)  {}

    /** */
    ~IterableFile () {}

    /** */
    dp::Iterator<Item>* iterator ()  { return new IteratorFile<Item> (_filename, _cacheItemsNb); }

    /** */
    int64_t getNbItems ()   {  return system::impl::System::file().getSize(_filename) / sizeof(Item); }

private:
    std::string     _filename;
    size_t          _cacheItemsNb;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_ */
