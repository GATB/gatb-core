/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
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
    IteratorFile (const std::string& filename, size_t cacheItemsNb=10000) :
        _filename(filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(cacheItemsNb), _isDone(true)

    {
        _file    = system::impl::System::file().newFile (filename, "rb");
        _buffer  = malloc (sizeof(Item) * _cacheItemsNb);
    }

    /** Destructor. */
    ~IteratorFile ()
    {
        if (_file)  { delete _file;  }
        if (_buffer) { free (_buffer); }
    }

    /** \copydoc Iterator::first */
    void first()
    {
        _cpt_buffer = 0;
        _idx        = 0;
        _isDone     = false;
        next ();
    }

    /** \copydoc Iterator::next */
    void next()
    {
        if (_cpt_buffer==0)
        {
            _idx = 0;
            _cpt_buffer = _file->fread (_buffer, sizeof(Item), _cacheItemsNb);
            if (_cpt_buffer==0)  { _isDone = true; }
        }

        this->_item =  ((Item*)_buffer) + _idx;
        _cpt_buffer --;
        _idx ++;
    }

    /** \copydoc Iterator::isDone */
    bool isDone()  { return _isDone; }

    /** \copydoc Iterator::item */
    Item& item ()  { return *(this->_item); }

    /** */
    size_t fill (std::vector<Item>& vec)
    {
        u_int64_t len = _file->getSize() / sizeof (Item);
        vec.resize(len);
        return _file->fread (vec.data(), sizeof(Item), len);
    }

private:
    std::string     _filename;
    system::IFile*  _file;
    void*           _buffer;
    int             _cpt_buffer;
    int             _idx;
    size_t          _cacheItemsNb;
    bool            _isDone;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_ */
