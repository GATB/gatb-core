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


/** \file Ukl_loader.hpp
 *  \date 03/04/2017
 *  \author grizk
 *  \brief Ukl_loader: universal kmer list loader
 */

#ifndef _GATB_CORE_KMER_IMPL_UKL_LOADER_HPP_
#define _GATB_CORE_KMER_IMPL_UKL_LOADER_HPP_

/********************************************************************************/

#include <iostream>


#include <stdio.h>
#include <climits>
#include <stdlib.h>
#include <math.h>

#include <assert.h>
#include <sys/time.h>
#include <string.h>
#include <fstream>



/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/



//simple bitArray with save/load from file
class simpleBitArray {
	
public:
	
	simpleBitArray() : _size(0)
	{
		_bitArray = nullptr;
	}
	
	simpleBitArray(uint64_t n) : _size(n)
	{
		_nchar  = (1ULL+n/64ULL);
		_bitArray =  (uint64_t *) calloc (_nchar,sizeof(uint64_t));
	}
	
	~simpleBitArray()
	{
		if(_bitArray != nullptr)
			free(_bitArray);
	}
	
	//copy constructor
	simpleBitArray(simpleBitArray const &r)
	{
		_size =  r._size;
		_nchar = r._nchar;
		_bitArray = (uint64_t *) calloc (_nchar,sizeof(uint64_t));
		memcpy(_bitArray, r._bitArray, _nchar*sizeof(uint64_t) );
	}
	
	void resize(uint64_t newsize)
	{
		_nchar  = (1ULL+newsize/64ULL);
		_bitArray = (uint64_t *) realloc(_bitArray,_nchar*sizeof(uint64_t));
		_size = newsize;
	}
	
	size_t size() const
	{
		return _size;
	}
	
	//clear whole array
	void clear()
	{
		memset(_bitArray,0,_nchar*sizeof(uint64_t));
	}
	
	//return value at pos
	uint64_t operator[](uint64_t pos) const
	{
		return (_bitArray[pos >> 6ULL] >> (pos & 63 ) ) & 1;
	}
	
	uint64_t get(uint64_t pos) const
	{
		return (*this)[pos];
	}
	
	//set bit pos to 1
	void set(uint64_t pos)
	{
		assert(pos<_size);
		_bitArray [pos >> 6] |=   (1ULL << (pos & 63) ) ;
		//__sync_fetch_and_or (_bitArray + (pos >> 6ULL), (1ULL << (pos & 63)) );
	}
	
	//set bit pos to 0
	void reset(uint64_t pos)
	{
		_bitArray [pos >> 6] &=   ~(1ULL << (pos & 63) ) ;
		//__sync_fetch_and_and (_bitArray + (pos >> 6ULL), ~(1ULL << (pos & 63) ));
	}
	
	void save(std::ostream& os) const
	{
		os.write(reinterpret_cast<char const*>(&_size), sizeof(_size));
		os.write(reinterpret_cast<char const*>(&_nchar), sizeof(_nchar));
		os.write(reinterpret_cast<char const*>(_bitArray), (std::streamsize)(sizeof(uint64_t) * _nchar));
	}
	
	void load(std::istream& is)
	{
		is.read(reinterpret_cast<char*>(&_size), sizeof(_size));
		is.read(reinterpret_cast<char*>(&_nchar), sizeof(_nchar));
		this->resize(_size);
		this->clear();
		is.read(reinterpret_cast<char *>(_bitArray), (std::streamsize)(sizeof(uint64_t) * _nchar));
	}
	
protected:
	uint64_t*  _bitArray;
	uint64_t _size;
	uint64_t _nchar;
};

	
	class UKL_loader
	{
		
	public:
		UKL_loader() : _ukl_set(0) {};
		~UKL_loader();
		
		simpleBitArray *  loadSetFromFile(int kmerSize, int mmerSize);
		simpleBitArray * getUKL();
		
	private:
		
		simpleBitArray * _ukl_set;
	};
	
	

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_UKL_LOADER_HPP_ */
