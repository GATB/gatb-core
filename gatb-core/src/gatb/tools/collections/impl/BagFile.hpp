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

/** \file BagFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bag implementation for file
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_FILE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>


#include <zlib.h>


/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/** Implementation of API from collections */
namespace impl          {
/********************************************************************************/

/** \brief Bag implementation for file
 */
template <typename Item> class BagFile : public Bag<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    BagFile (const std::string& filename) : _filename(filename), _file(0)
    {
        /** We first erase the file. */
        system::impl::System::file().remove (filename);

        /** We get a handle on the file. */
        _file = system::impl::System::file().newFile (filename, "wb+");
    }

    /** Destructor. */
    ~BagFile ()
    {
        if (_file)  { delete _file; }
    }

    /** */
    const std::string& getName () const { return _filename; }

    /**  \copydoc Bag::insert */
    void insert (const Item& item)  {  _file->fwrite (&item, sizeof(Item), 1);  }

    void insert (const std::vector<Item>& items, size_t length)
    {
        if (length == 0)  { length = items.size(); }
        _file->fwrite (items.data(), sizeof(Item), length);
    }

    /** */
    void insert (const Item* items, size_t length)
    {
        _file->fwrite (items, sizeof(Item), length);
    }

    /**  \copydoc Bag::flush */
    void flush ()  { _file->flush(); }

private:
    std::string _filename;
    system::IFile* _file;
};

    

/** \brief Bag implementation for file
 */
template <typename Item> class BagGzFile : public Bag<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    BagGzFile (const std::string& filename) : _filename(filename), _gzfile(0)
    {
        /** We first erase the file. */
        system::impl::System::file().remove (filename);
        
        /** We get a handle on the file. */
        _gzfile =   gzopen(filename.c_str(),"wb1"); // (gzFile *)
    }
    
    /** Destructor. */
    ~BagGzFile ()
    {
        if (_gzfile)  {  gzclose(_gzfile); }
    }
    
    /** */
    const std::string& getName () const { return _filename; }
    
    /**  \copydoc Bag::insert */
    void insert (const Item& item)  { gzwrite(_gzfile,&item,sizeof(Item)); }
    
    void insert (const std::vector<Item>& items, size_t length)
    {
        if (length == 0)  { length = items.size(); }
        //_file->fwrite (items.data(), sizeof(Item), length);
        gzwrite(_gzfile,items.data(),sizeof(Item)*length);
    }
    
    /** */
    void insert (const Item* items, size_t length)
    {
        gzwrite(_gzfile,items,sizeof(Item)*length);
    }
    
    /**  \copydoc Bag::flush */
    void flush ()  {
        gzflush(_gzfile,Z_FINISH); // not sure if necessary
    }
    
private:
    std::string _filename;
    system::IFile* _file;
    gzFile  _gzfile;
};

    
//stores items in a file, compressed by storing count, and only diff with previous item
// todo writing is not yet cached
    //todo implem codage diff
template <typename Item> class BagCountCompressedFile : public Bag<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    BagCountCompressedFile (const std::string& filename) : _filename(filename), _file(0),_previous()
    {
        /** We first erase the file. */
        system::impl::System::file().remove (filename);
        _buffer_size = 16384;
        /** We get a handle on the file. */
        _file = system::impl::System::file().newFile (filename, "wb+");
        _bufferOut  = (u_int8_t*) malloc (sizeof(u_int8_t) * _buffer_size);
        _idx=0;
        _size_item = sizeof(Item);

    }
    
    /** Destructor. */
    ~BagCountCompressedFile ()
    {
        if (_file)  { delete _file; }
        free( _bufferOut );
    }
    
    /** */
    const std::string& getName () const { return _filename; }
    
    /**  \copydoc Bag::insert */
    void insert (const Item& item)  {
        printf("insert single element not supported in BagCountCompressedFile\n");

    }
    
    void insert (const std::vector<Item>& items, size_t length)
    {
        if (length == 0)  { length = items.size(); }
        
        this->insert(items.data(), length);
        
    }
    
    /** */
    void insert (const Item* items, size_t length)
    {
      //  Item diff_to_prev ;
        _sizeOutput =0;
        
        _previous = items[0];
        u_int8_t abundance = 0;

        for (int ii=0; ii< length; ii++) {
            //diff_to_prev =  items[i] - _previous;
            if (items[ii] == _previous)  {
                
                abundance++;

                if(abundance == UINT8_MAX)
                {
                    write(abundance,_previous);
                    abundance     = 0;
                }
                
            }
            else
            {
                write(abundance,_previous);
                abundance     = 1;
                _previous = items[ii];
            }
        }
        //last item
        write(abundance,_previous);


        
        //printf("In %lu B  Out %llu  B  ratio  %f \n",length*sizeof(Item), _sizeOutput, length*sizeof(Item) / (float) _sizeOutput);
    }
    
    /**  \copydoc Bag::flush */
    void flush ()  {
        //printf("--- flush  BagCountCompressedFile %d this %p ---\n",_idx,this);
        flushBuffer ();
        _file->flush();
    }
    
    //maybe use a bagfile encapsulated in a bagcache to write big chunks to file..
    //or internal small buffer  16KB
private:
    void write (u_int8_t  abundance,Item  elem)  {
        // printf("writing abundance %i \n",abundance);
        
        if( (_idx + (1+_size_item)) >= _buffer_size )
        {
            flushBuffer ();
        }
        memcpy(_bufferOut+_idx,&abundance, 1 ); _idx+= 1;
        memcpy(_bufferOut+_idx,&elem, _size_item ); _idx+= _size_item;
        
//        _file->fwrite (&abundance, sizeof(u_int8_t), 1);
//        _file->fwrite (&elem, sizeof(Item), 1);
        _sizeOutput += 1+ sizeof(Item);
    }

    void flushBuffer ()  {
        _file->fwrite (_bufferOut, sizeof(u_int8_t), _idx);
        _idx =0;

    }
    
    std::string _filename;
    system::IFile* _file;
    Item _previous;
    u_int64_t _sizeOutput;
    u_int8_t*           _bufferOut;
    int             _idx;
    size_t _size_item;
    size_t _buffer_size;

};
    
    
/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_FILE_HPP_ */
