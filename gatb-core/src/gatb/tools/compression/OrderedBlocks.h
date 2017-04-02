/*****************************************************************************
 *   Leon: reference free compression for NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Benoit, G.Rizk, C.Lemaitre
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

#ifndef __bloocoo_ReadWriter__
#define __bloocoo_ReadWriter__

/********************************************************************************/

#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/bank/impl/Banks.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include <pthread.h>


/********************************************************************************/
/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::bank::impl;
using namespace gatb::core::tools::collections;

/********************************************************************************/

class  OrderedBlocks;

typedef struct
{
     OrderedBlocks * obw;
} thread_args;

void * writer(void * args);





class  OrderedBlocks : public SmartPointer
{
public:
    
    /** Constructor. */
     OrderedBlocks (IFile*  outputfile, size_t buffsize)
    :   _writer_available(1),  _buffer_full(0),_to_be_written(0) ,_buffWrite (buffsize), _buffReceive (buffsize),_ref(0),_os(0)
	, _nbMax(buffsize), _base(0), _idx(0)
    {
        setRef(outputfile);
        
        //create writer thread
        t_arg.obw = this;
		
        pthread_mutex_init(&writer_mutex, NULL);
        pthread_cond_init (&writer_available_cond, NULL);
        pthread_cond_init (&buffer_full_cond, NULL);
        
        pthread_create (&_thread, NULL,  writer, &t_arg);
		//  printf("  .. end   this %p    cond %p   targs %p\n",this,&buffer_full_cond, &t_arg );
		
    }
	
	
	//constructor with output stream
	OrderedBlocks (tools::storage::impl::Storage::ostream * outstream, size_t buffsize)
	:   _writer_available(1),  _buffer_full(0),_to_be_written(0) ,_buffWrite (buffsize), _buffReceive (buffsize),_ref(0),_os(0)
	, _nbMax(buffsize), _base(0), _idx(0)
	{
		_os = outstream;
		t_arg.obw = this;
		pthread_mutex_init(&writer_mutex, NULL);
		pthread_cond_init (&writer_available_cond, NULL);
		pthread_cond_init (&buffer_full_cond, NULL);
		pthread_create (&_thread, NULL,  writer, &t_arg);
		
	}
	
    /** Destructor. */
    ~ OrderedBlocks ()
    {
        setRef(0);
    }
    
    
    void insert (u_int8_t* data, u_int64_t size, int blockId) ;
    
    void incDone (int nbdone);
	
    void waitForWriter ();
    void FlushWriter ();
	
    
    thread_args t_arg;
    pthread_mutex_t writer_mutex;
    pthread_cond_t writer_available_cond;
    pthread_cond_t buffer_full_cond;
    int _writer_available ;
    int _buffer_full;
    int _to_be_written;
    
	
	std::vector<   std::vector<u_int8_t>    >      _buffWrite;
	std::vector<  std::vector<u_int8_t>    >      _buffReceive;
	
	
    
	
	IFile* _ref;
	
	gatb::core::tools::storage::impl::Storage::ostream * _os;
	
	

	
protected:
    
    
    void setRef (IFile* ref)
    {
        _ref=ref;
    }
	
	
	
    //they will be swapped _buffWrite.swap(_buffReceive);
    
	
	
	
    size_t      _nbMax; // buffer size
	size_t      _base;
    size_t      _idx;   //number of sequences stored in the buffer
	
    pthread_t  _thread; // the writer thread
	
	
    
private:
    //assign
     OrderedBlocks& operator=(const  OrderedBlocks& bk){
        return *this;
    }
    //copy construct
     OrderedBlocks(const  OrderedBlocks& bk){ }
	
    
};




#endif /* defined(__bloocoo_xcode__ReadWriter__) */
