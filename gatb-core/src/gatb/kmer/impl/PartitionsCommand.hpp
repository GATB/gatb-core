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

#ifndef _PARTITIONSCOMMAND__HPP_
#define _PARTITIONSCOMMAND__HPP_

/********************************************************************************/

#define IX(x,rad) ((rad)+(256)*(x))


#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/bank/impl/Banks.hpp>

#include <queue>
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::kmer::impl;


/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
	
	class MemAllocator;

	//class containing info of each parti : exact number of kmers per parti
	//will be computed by fillparti, then used by fillsolids
	//needed beacause the exact number of kmers per parti can no longer be inferred from partition file size
	// xmer ==1 -->  only k0 mer    (xmer == nb max of kmers in kx mer, not the x   k1 mer == 2 kmers max in k1 mer)
	template <size_t xmer>
	class PartiInfo
	{
	public:
		
		inline void incKmer(int numpart, u_int64_t val=1)
		{
			_nb_kmers_per_parti[numpart]+=val;
		}
		
		inline void incKxmer(int numpart, u_int64_t val=1)
		{
			_nb_kxmers_per_parti[numpart]+=val; //now used to store number of kxmers per parti
		}
		
		//superksize in number of kmer inside this superk
		//increases both superk and kmer count
		inline void incSuperKmer_per_minimBin(int numbin, int superksize, u_int64_t val=1)
		{
			_superk_per_mmer_bin[numbin]+= val ;
			_kmer_per_mmer_bin[numbin]+= val *superksize;
			
		}
		
		//kxmer count (regardless of x size), used for ram
		inline void incKxmer_per_minimBin(int numbin, u_int64_t val=1)
		{
			_kxmer_per_mmer_bin[numbin]+= val ;
		}
		
		//numaprt, radix, size of kx mer
		inline void incKmer_and_rad(int numpart, int radix,int x,  u_int64_t val=1) //appele ds vrai loop
		{
			_nb_kxmers_per_parti[numpart] += val; //now used to store number of kxmers per parti
			_nb_kmers_per_parti[numpart]+= (val * (x+1)); // number of  'real' kmers
			_nbk_per_radix_per_part[x][radix][numpart]+=val; // contains number of kx mer per part per radix per x
		}
		
		
		PartiInfo& operator+=  (const PartiInfo& other)
		{
			
			//add other parti info , synced
			
			for (int np=0; np<_nbpart; np++) {
				
				for (int xx=0; xx<xmer; xx++)
				for (int rad=0; rad<256; rad++) {
					__sync_fetch_and_add( & (_nbk_per_radix_per_part[xx][rad][np]),   other.getNbKmer(np,rad,xx) );
				}
			
				__sync_fetch_and_add( _nb_kmers_per_parti + np,        other.getNbKmer(np) );
				__sync_fetch_and_add( _nb_kxmers_per_parti + np,   other.getNbSuperKmer(np) );
				
			}
			
			for (int ii=0; ii< _num_mm_bins; ii++) {
				__sync_fetch_and_add( _superk_per_mmer_bin + ii,   other.getNbSuperKmer_per_minim(ii) );
				__sync_fetch_and_add( _kmer_per_mmer_bin + ii,   other.getNbKmer_per_minim(ii) );
				__sync_fetch_and_add( _kxmer_per_mmer_bin + ii,   other.getNbKxmer_per_minim(ii) );
				
			}
			
			return *this;
		}
		
		
		inline  u_int64_t getNbKmer(int numpart) const
		{
			return _nb_kmers_per_parti[numpart];
		}
		
		//get nbk in bin radix of parti numpart
		inline  u_int64_t getNbKmer(int numpart, int radix, int xx) const
		{
			return _nbk_per_radix_per_part[xx][radix][numpart];
		}
		
		
		inline  u_int64_t   getNbSuperKmer(int numpart) const //now used for number of kxmers (indistinctive of xsize)
		{
			return _nb_kxmers_per_parti[numpart];
		}
		
		
		inline  u_int64_t   getNbSuperKmer_per_minim(int numbin) const
		{
			return _superk_per_mmer_bin[numbin];
		}
		
		inline  u_int64_t   getNbKmer_per_minim(int numbin) const
		{
			return _kmer_per_mmer_bin[numbin];
		}
		
		inline  u_int64_t   getNbKxmer_per_minim(int numbin) const
		{
			return _kxmer_per_mmer_bin[numbin];
		}
		
		
		void clear()
		{
			memset(_nb_kmers_per_parti, 0, _nbpart*sizeof(u_int64_t));
			memset(_nb_kxmers_per_parti, 0, _nbpart*sizeof(u_int64_t));
			memset(_superk_per_mmer_bin, 0, _num_mm_bins *sizeof(u_int64_t));
			memset(_kmer_per_mmer_bin, 0, _num_mm_bins *sizeof(u_int64_t));
			memset(_kxmer_per_mmer_bin, 0, _num_mm_bins *sizeof(u_int64_t));

			
			for (int xx=0; xx<xmer; xx++)
			for(int ii=0; ii<256; ii++)
			{
				memset(_nbk_per_radix_per_part[xx][ii], 0, _nbpart*sizeof(u_int64_t));
			}
			
		}
		
		void printInfo() const
		{
			
			printf("------------------\n");
			printf("Nb kmers per parti\n");
			
			for (int np=0; np<_nbpart; np++) {
				printf("Parti[%i]= %lli\n",np,this->getNbKmer(np));
			}
			
			
			printf("------------------------\n");
			printf("Nb kxmers per parti\n");
			
			for (int np=0; np<_nbpart; np++) {
				printf("Parti[%i]= %lli\n",np,this->getNbSuperKmer(np));
			}
			
			
//			printf("----------------------------\n");
//			printf("Nb kmers per parti per radix\n");
//			
//			for (int np=0; np<_nbpart; np++) {
//				printf("___ Parti %i ___\n",np);
//				
//				for (int rad=0; rad<256; rad++) {
//					printf("%10lli  ",this->getNbKmer(np,rad,0));
//					if((rad & 7) == 7) printf("\n");
//				}
//				printf("\n");
//				
//			}
			
			printf("----------------------------\n");
			printf("Nb Super kmers , nb kmers per minim bin\n");
			
			u_int64_t sumk = 0;
			u_int64_t sumsuperk = 0;
			
			for (int np=0; np<_num_mm_bins; np++) {
				typedef typename Kmer<31>::Type           Typem; //should be kmer size 
				Typem cur = np;
				
				printf("Bin[%5i (%s) ]= %lli    %lli\n",np,cur.toString(_mm).c_str(), this->getNbSuperKmer_per_minim(np),this->getNbKmer_per_minim(np)    );
				
				sumk += this->getNbKmer_per_minim(np);
				sumsuperk +=  this->getNbSuperKmer_per_minim(np);
			}
			
			printf("total number of kmers %lli  total number of superkmers %lli \n",sumk,sumsuperk);
			printf("Average size of superkmers :  %f \n",(double)sumk / (double)sumsuperk );
			
			
			
		}
		
		
		PartiInfo(int nbpart, int minimsize) : _nbpart(nbpart), _mm(minimsize)
		{
			_nb_kmers_per_parti =  (u_int64_t  *)  calloc(nbpart,sizeof(u_int64_t));
			_nb_kxmers_per_parti =  (u_int64_t  *)  calloc(nbpart,sizeof(u_int64_t));
			
			_num_mm_bins =   1 << (2*_mm);
			_superk_per_mmer_bin = ( u_int64_t * )  calloc (_num_mm_bins ,sizeof(u_int64_t) );
			_kmer_per_mmer_bin = ( u_int64_t * )  calloc (_num_mm_bins ,sizeof(u_int64_t) );
			_kxmer_per_mmer_bin = ( u_int64_t * )  calloc (_num_mm_bins ,sizeof(u_int64_t) );


			
			for(int xx=0; xx<xmer; xx++)
			for(int ii=0; ii<256; ii++)
			{
				_nbk_per_radix_per_part[xx][ii] = (u_int64_t  *) calloc(nbpart,sizeof(u_int64_t));
			}
			
		}
		
		
		PartiInfo(const PartiInfo& cr) //need copy constr,  (for fillparti class par exemple)
		//the copy contr realloc its own  arrays, zero init
		{
			
			_num_mm_bins = cr._num_mm_bins;
			_nbpart =  cr._nbpart;
			_mm = cr._mm;
			
			
			_nb_kmers_per_parti =  (u_int64_t  *)  calloc(_nbpart,sizeof(u_int64_t));
			_nb_kxmers_per_parti =  (u_int64_t  *)  calloc(_nbpart,sizeof(u_int64_t));
			_superk_per_mmer_bin = ( u_int64_t * )  calloc (_num_mm_bins ,sizeof(u_int64_t) );
			_kmer_per_mmer_bin = ( u_int64_t * )  calloc (_num_mm_bins ,sizeof(u_int64_t) );
			_kxmer_per_mmer_bin = ( u_int64_t * )  calloc (_num_mm_bins ,sizeof(u_int64_t) );
			
			for(int xx=0; xx<xmer; xx++)
			for(int ii=0; ii<256; ii++)
			{
				_nbk_per_radix_per_part[xx][ii] = (u_int64_t  *) calloc(_nbpart,sizeof(u_int64_t));
			}
			
			//	printf("PartiInfo copy constr %p  _nb_kmers_per_parti %p \n",this,_nb_kmers_per_parti);
			
		}
		
		~PartiInfo()
		{
			free(_nb_kmers_per_parti);
			free(_nb_kxmers_per_parti);
			free(_superk_per_mmer_bin);
			free(_kmer_per_mmer_bin);
			free(_kxmer_per_mmer_bin);
			
			for(int xx=0; xx<xmer; xx++)
			for(int ii=0; ii<256; ii++)
			{
				free(_nbk_per_radix_per_part[xx][ii]);
			}
			
			//	printf("print info destroyed %p  _nb_kmers_per_parti %p \n",this,_nb_kmers_per_parti);
		}
		
	private:
		
		u_int64_t  * _nb_kmers_per_parti;
		u_int64_t  * _nb_kxmers_per_parti; //now used to store number of kxmers per parti
		u_int64_t * _superk_per_mmer_bin;
		u_int64_t * _kmer_per_mmer_bin;
		u_int64_t * _kxmer_per_mmer_bin;

		u_int64_t * _nbk_per_radix_per_part[xmer][256];//number of kxmer per parti per rad
		u_int64_t _num_mm_bins;
		
		int _nbpart;
		int _mm;
		
	};
	
	
	

/********************************************************************************/

	
/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
class PartitionsCommand : public ICommand, public system::SmartPointer
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;

	

	
	
	
    PartitionsCommand (
        Bag<Count>* solidKmers,
        Iterable<Type>&   partition,
        IHistogram*    histogram,
        ISynchronizer* synchro,
        u_int64_t&     totalKmerNbRef,
        size_t         abundance,
        IteratorListener* progress,
		PartiInfo<5> * pInfo,
		int parti,
		size_t      nbCores,
		size_t      kmerSize,
	    MemAllocator * pool


    );

    ~PartitionsCommand();

protected:
    size_t              _abundance;
    BagCache<Count>      _solidKmers;
    Iterable<Type>&      _partition;
    HistogramCache      _histogram;
    ProgressSynchro     _progress;
    u_int64_t           _totalKmerNb;
    u_int64_t&          _totalKmerNbRef;
	PartiInfo<5> * _pInfo;
	int _parti_num;
    size_t      _nbCores;
	size_t      _kmerSize;
	MemAllocator * _pool;
	
    void insert (const Count& kmer);
};

/********************************************************************************/
/** */
template<size_t span>
class PartitionsByHashCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;


    PartitionsByHashCommand (
        Bag<Count>*      solidKmers,
        Iterable<Type>&  partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        size_t          abundance,
        IteratorListener* progress,
        u_int64_t       hashMemory,
		PartiInfo<5> * pInfo,
		int parti,
		size_t      nbCores,
		size_t      kmerSize,
	    MemAllocator * pool

    );

    void execute ();

private:
    u_int64_t _hashMemory;
};

		
/********************************************************************************/
/** */
template<size_t span>
class PartitionsByVectorCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;


	//used for the priority queue
	typedef std::pair<int, Type> kxp; //id pointer in vec_pointer , value
	
	struct kxpcomp {
		bool operator() (kxp l,kxp r) { return ((r.second) < (l.second)); }
	} ;
	

    PartitionsByVectorCommand (
        Bag<Count>*  solidKmers,
        Iterable<Type>&    partition,
        IHistogram*     histogram,
        ISynchronizer*  synchro,
        u_int64_t&      totalKmerNbRef,
        size_t          abundance,
        IteratorListener* progress,
		PartiInfo<5> * pInfo, //x+1  (space for k+0,  .. ,k+x mer)
		int parti,
		size_t      nbCores,
		size_t      kmerSize,
	    MemAllocator * pool
    );

    void execute ();

private:

//    vector<Type> kmers;
	
//	vector <   vector < vector<Type> >  >    radix_kmers; //by  xmer by radix 4nt bins //pas tres beau ,bcp indirection
	Type ** _radix_kmers;

};

	
template<size_t span>
class SortCommand : public ICommand, public system::SmartPointer
{
public:
	typedef typename Kmer<span>::Type  Type;

	SortCommand( Type **   kmervec, int begin, int end, uint64_t * radix_sizes) : _radix_kmers(kmervec),_deb(begin), _fin(end), _radix_sizes(radix_sizes) {}
	
	 void execute ()
	{

		//printf("will exec  range [ %i ; %i] \n ",_deb,_fin);
		for (int ii= _deb; ii<= _fin; ii++) {
	
//				if(_radix_kmers[ii].size() > 0)
//					std::sort (_radix_kmers[ii].begin (), _radix_kmers[ii].end ());
			
			if(_radix_sizes[ii] > 0)
				std::sort (&_radix_kmers[ii][0] , &_radix_kmers[ii][ _radix_sizes[ii]]);
		}
		
	}
	
	private :
	
	int _deb, _fin;
	//std::vector <  vector<Type> >   & _radix_kmers;
	Type ** _radix_kmers;
	uint64_t * _radix_sizes;
};
	

	

template<size_t span>
class KxmerPointer
{
public:
	typedef typename Kmer<span>::Type  Type;
	
	//on lui passe le vector dun kxmer //std::vector<vector<Type> >  & kmervec
	KxmerPointer(Type ** kmervec, int prefix_size, int x_size, int min_radix, int max_radix, int kmerSize, uint64_t * radix_sizes) : _kxmers(kmervec)
	,_prefix_size(prefix_size), _x_size(x_size),_cur_idx(-1) ,_kmerSize(kmerSize),_low_radix(min_radix),_high_radix(max_radix), _radix_sizes(radix_sizes)  {
		
		_idx_radix = min_radix;
		Type un = 1;
		_kmerMask = (un << (_kmerSize*2)) - un;
		
		_shift_size = ( (4 - _prefix_size) *2) ;
		_radixMask = Type(_idx_radix) ;
		_radixMask = _radixMask << ((_kmerSize-4)*2);
		_radixMask = _radixMask  << (2*_prefix_size)  ;

	}
	
	
	
	inline bool next ()
	{

		_cur_idx++;
		
		//go to next non empty radix
	//	while(_idx_radix<= _high_radix &&_cur_idx >= _kxmers[_idx_radix].size()  ) //if using vectors
		while(_idx_radix<= _high_radix &&_cur_idx >=   _radix_sizes[_idx_radix]    )
		{
			_idx_radix++;
			_cur_idx = 0;
			//update radix mask does not happen often
			_radixMask = Type(_idx_radix) ;
			_radixMask = _radixMask << ((_kmerSize-4)*2);
			_radixMask = _radixMask  << (2*_prefix_size)  ;
		}
		
		return (_idx_radix <= _high_radix);
	}
	
	inline Type value() const
	{
		Type res =  ( ((_kxmers[_idx_radix][_cur_idx]) >> _shift_size)  |  _radixMask  ) & _kmerMask ;

		return res ;
	}
	
	private :
	
	//std::vector <  vector < Type > >    & _kxmers;
	Type **    _kxmers;
	uint64_t * _radix_sizes;

	int64_t _cur_idx;
	Type _cur_radix;
	Type _kmerMask;
	Type _radixMask;
	
	int _idx_radix;
	int _low_radix, _high_radix;
	int _shift_size;
	int _prefix_size;
	int _kmerSize;
	int _x_size; //x size of the _kxmersarray
	
	
};


	
//make it an allocator usable by std vector ?
class MemAllocator
{
	public:
	
	
	//clear all previous allocs, and alloc pool capacity
	void reserve(u_int64_t size)
	{
		//printf("request pool reserve %llu B \n",size);
		if(size ==0 && mainbuffer !=NULL)
		{
			free(mainbuffer);
			capacity = used_space = 0;
			mainbuffer = NULL ;
		}
		
		mainbuffer =(char *)  malloc(size);
		capacity  = size;
		used_space = 0;
	}
	
	
	//should be thread safe
	char * pool_malloc(u_int64_t requested_size)
	{
		
		u_int64_t synced_used_space = __sync_fetch_and_add(&used_space, requested_size);

		if(requested_size> (capacity - synced_used_space))
		{
			printf("mem allocator full\n");
			__sync_fetch_and_add(&used_space, -requested_size);
			return NULL;
		}
		
		
		return mainbuffer + synced_used_space;
	}
	
	u_int64_t getCapacity()
	{
		return capacity;
	}
	
	u_int64_t getUsedSpace()
	{
		return used_space;
	}

	
	void free_all()
	{
		used_space = 0;
	}
	
	MemAllocator() : capacity(0),used_space(0),mainbuffer(NULL)
	{
		
	}

	
	~MemAllocator()
	{
		if(mainbuffer !=NULL)
			free(mainbuffer);
	}
		
	private :
	
	
	char * mainbuffer;
	u_int64_t capacity; //in bytes
	u_int64_t used_space;
	
		
};
	
	

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/


#endif
