
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


#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>
#include <gatb/tools/math/Integer.hpp>


#include "OrderedBlocks.h"

#define DEBUG(a)  //printf a


// We use the required packages
using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;
using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;
using namespace gatb::core::tools::math;

void  OrderedBlocks::insert(u_int8_t* data, u_int64_t size, int blockId)
{
	DEBUG(("___ oblock insert bid %i \n",blockId));

    size_t id ;
	
    while (1)
    {
		id  = blockId - _base;
        if(id < _nbMax )
        {
			//cout << blockId << " " << id << " " << size << endl;
			_buffReceive[id] =  std::vector<u_int8_t> (data,data +  size);
			DEBUG(("___ successfully  inserted bid %i  rid %i  _nbMax %i  _idx %i   base %i \n",blockId,id,_nbMax,_idx,_base));

            break;
        }
        else
        {
            //the current thread has finished the block and is now ahead of others,
            //wait for other threads to finish filling the buffer
            
            pthread_mutex_lock(&writer_mutex);
			id  = blockId - _base; // id value must be update here, because can be changed in between
            while (id >= _nbMax )
            {
				DEBUG(("___________ waiting to insert   bid %i   id %i, _nbMax %i  _idx %i   _base %i   going to sleep\n",blockId,id,_nbMax,_idx,_base));

                pthread_cond_wait(&buffer_full_cond, &writer_mutex);
				id  = blockId - _base;
				DEBUG(("___ waken up  to insert   bid %i   id %i, _nbMax %i  _base %i\n",blockId,id,_nbMax,_base));

            }
            pthread_mutex_unlock(&writer_mutex);
            
			
        }
    }
    
    
}


void  OrderedBlocks::incDone (int nbdone)
{
    
    size_t old_val;
	if(nbdone==0)
		return;
	
    old_val = __sync_fetch_and_add (& _idx, nbdone);
	 DEBUG(("\n ++ ninc oldval  %i base = %i ++ \n",old_val,_base));
	
    // buffer completely filled : activates writer thread
    if( (old_val + nbdone) ==_nbMax)
    {
		 DEBUG(("\n ++ ninc done %zd  and buff done, base = %zd ++ \n",old_val,_base));
		
        ///// wait for writer to be free
        pthread_mutex_lock(&writer_mutex);
        while (_writer_available==0) {
			DEBUG(("\n ++ wait for writer to be free ++ \n"));

            pthread_cond_wait(&writer_available_cond, &writer_mutex);
        }
        //buffReceive is full, will be written
		
		//cout << "lala: " << _buffReceive.size() << endl;

        _buffWrite.swap(_buffReceive);
        
        _to_be_written = _nbMax;
        _idx=0;
        _base += _nbMax;

        //signal writer he should write buffer
        _buffer_full = 1;
		DEBUG(("\n ++ swap buffers , order writer to write++ idx %i new base %i should wake all  \n",_idx,_base));

        pthread_cond_broadcast(&buffer_full_cond);
        pthread_mutex_unlock(&writer_mutex);
		
		
        
    }
    
}

void  OrderedBlocks::waitForWriter ()
{
    pthread_mutex_lock(&writer_mutex);
    while (_writer_available==0) {
		//     printf("worker going to sleep, waiting for writer to be finished .. \n");
        pthread_cond_wait(&writer_available_cond, &writer_mutex);
    }
    pthread_mutex_unlock(&writer_mutex);
	//   printf("end waitForWriter \n");
	
}




//guaranteed by design only one thread at the end will call this
void  OrderedBlocks::FlushWriter ()
{
	pthread_mutex_lock(&writer_mutex);

    if( _idx)
    {
		  DEBUG(("\n flushing %zd base = %zd ++ \n",_idx,_base));
        ///// wait for writer to be free
        while (_writer_available==0) {
            pthread_cond_wait(&writer_available_cond, &writer_mutex);
        }
		
        
        _buffWrite.swap(_buffReceive);
        _to_be_written = _idx;
        
        //signal writer he should write buffer
        _buffer_full = 1;
        pthread_cond_broadcast(&buffer_full_cond);
		
    }
    
	pthread_mutex_unlock(&writer_mutex);

	
    //wait again for writer to finish flushing
    pthread_mutex_lock(&writer_mutex);
    while (_buffer_full==1) {
        pthread_cond_wait(&writer_available_cond, &writer_mutex);
    }
    pthread_mutex_unlock(&writer_mutex);
    
    //cout << "FLUSHED" << endl;
    _base = 0;
    _writer_available = 1;
    _buffer_full = 0;
    _idx = 0;
    _to_be_written = 0;
    //_nbMax(buffsize), _buffWrite (buffsize), _buffReceive (buffsize), _base(0), _writer_available(1), _buffer_full(0), _idx(0),_to_be_written(0)
}





// writer thread gets a pointer to  orderbankwriter object
void * writer(void * args)
{
	
    thread_args *targ = (thread_args*) args;
     OrderedBlocks * obw = targ->obw;
	
    IFile* outbank = obw->_ref;
	
	tools::storage::impl::Storage::ostream * outstream =  obw->_os;

	
    int  * to_be_written =  &(obw->_to_be_written);
    std::vector< std::vector<u_int8_t> > * _buffWrite = &(obw->_buffWrite);
	
    
	
    //when waken, writes content of  vector  _buffWrite to bank
    while(1)
    {
        pthread_mutex_lock(&(obw->writer_mutex));
        while (obw->_buffer_full==0)
        {
            pthread_cond_wait(&(obw->buffer_full_cond), &(obw->writer_mutex));
        }
        obw->_writer_available = 0; //writer becomes busy
		  DEBUG((" **** writer thread  awaken !! ..  will write %i elems **** \n",*to_be_written));
        pthread_mutex_unlock(&(obw->writer_mutex));
        
		
        //writes the buffer

		
        for ( std::vector< std::vector<u_int8_t>  >::iterator it = _buffWrite->begin(); (it != _buffWrite->end()) && (*to_be_written) ; it++)
        {
            //outbank->insert((*it));
            //cout << "Size: " << (*it).size() << endl;
			if(outstream==0 && outbank!=0)
				outbank->fwrite( &(*it)[0] ,(*it).size(), 1);

			// replace file by write to storage group/dataset through outstream
			if(outstream!=0 )
				outstream->write(reinterpret_cast<char const*>(&(*it)[0]), (*it).size());
			
			
            (*to_be_written) --;
        }
        
        
        //signal to others buffer has been written and writer is available again
        pthread_mutex_lock(&(obw->writer_mutex));
        obw->_buffer_full=0;
        obw->_writer_available=1;
		//  printf(" writer thread  setting  _buffer_full to 0 .. %i  tobewrit %i \n",obw->_buffer_full,*to_be_written);
        pthread_cond_broadcast(&(obw->writer_available_cond));
		DEBUG((" **** writer thread  finished !!            **** \n"));

        pthread_mutex_unlock(&(obw->writer_mutex));
    }
    
    
}

