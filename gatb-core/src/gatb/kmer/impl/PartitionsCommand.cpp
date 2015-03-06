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

#include <gatb/kmer/impl/PartitionsCommand.hpp>

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

#define DEBUG(a)  // printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

#define IX(x,rad) ((rad)+(256)*(x))

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
PartitionsCommand<span>:: PartitionsCommand (
    Bag<Count>*         solidKmers,
    Iterable<Type>&     partition,
    IHistogram*         histogram,
    ISynchronizer*      synchro,
    u_int64_t&          totalKmerNbRef,
    std::pair<size_t,size_t> abundance,
    IteratorListener*   progress,
    TimeInfo&           timeInfo,
    PartiInfo<5>&       pInfo,
    int                 parti,
    size_t              nbCores,
    size_t              kmerSize,
    MemAllocator&       pool,
    size_t              cacheSize
)
    : _abundance(abundance),
      _solidKmers(solidKmers, cacheSize, synchro),
      _partition(partition),
      _histogram(histogram),
      _progress(progress),
      _globalTimeInfo(timeInfo),
      _totalKmerNb(0),
      _totalKmerNbRef(totalKmerNbRef),
      _pInfo(pInfo),
      _parti_num(parti),
      _nbCores(nbCores),
      _kmerSize(kmerSize),
      _pool(pool)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
PartitionsCommand<span>::~PartitionsCommand()
{
    __sync_fetch_and_add (&_totalKmerNbRef, _totalKmerNb);

    _globalTimeInfo += _timeInfo;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsCommand<span>::insert (const Count& kmer)
{
    _totalKmerNb++;

    /** We should update the abundance histogram*/
    _histogram.inc (kmer.abundance);

    /** We check that the current abundance is in the correct range. */
    if (kmer.abundance >= this->_abundance.first && kmer.abundance <= this->_abundance.second)  {  this->_solidKmers.insert (kmer);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsCommand<span>::insert (const Type& kmer, const SolidityCounter& counter)
{
    _totalKmerNb++;

    /** Shortcut. */
    SolidityCounter::Int actualCount = counter.computeSum();

    /** We should update the abundance histogram*/
    _histogram.inc (actualCount);

    /** We check that the current abundance is in the correct range. */
    if (counter.isSolid () == true)  {  this->_solidKmers.insert (Count(kmer,actualCount));  }

    //if (actualCount >= this->_abundance && actualCount <= max_couv)  {  this->_solidKmers.insert (Count(kmer,actualCount));  }
}

/*********************************************************************
                #     #     #      #####   #     #
                #     #    # #    #     #  #     #
                #     #   #   #   #        #     #
                #######  #     #   #####   #######
                #     #  #######        #  #     #
                #     #  #     #  #     #  #     #
                #     #  #     #   #####   #     #
*********************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
/** in this scheme we count k-mers inside a partition by a hash table */
template<size_t span>
PartitionsByHashCommand<span>:: PartitionsByHashCommand (
    Bag<Count>*             solidKmers,
    Iterable<Type>&         partition,
    IHistogram*             histogram,
    ISynchronizer*          synchro,
    u_int64_t&              totalKmerNbRef,
    std::pair<size_t,size_t> abundance,
    IteratorListener*       progress,
    TimeInfo&               timeInfo,
    PartiInfo<5>&           pInfo,
    int                     parti,
    size_t                  nbCores,
    size_t                  kmerSize,
    MemAllocator&           pool,
    size_t                  cacheSize,
    u_int64_t               hashMemory
)
    : PartitionsCommand<span> (
        solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress, timeInfo, pInfo,parti,nbCores,kmerSize,pool,cacheSize),
        _hashMemory(hashMemory)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsByHashCommand<span>:: execute ()
{
	size_t count=0;

	/** We need a map for storing part of solid kmers. */
	OAHash<Type> hash (_hashMemory); //or use hash16 to ensure always finishes ?

	/** We directly fill the vector from the current partition file. */
	Iterator<Type>* it = this->_partition.iterator();  LOCAL(it);

	// If the partition holds kmers (and not superkmers), it would be :
	//      for (it->first(); !it->isDone(); it->next())   {  hash.increment (it->item());  }

	DEBUG (("PartitionsByHashCommand::execute:  fillsolid parti num %i  by oahash --- mem %llu  MB\n",
        this->_parti_num,_hashMemory/MBYTE
    ));
	
	//with decompactage
	//superk
	u_int8_t		nbK, rem ;
	Type compactedK;
	int ks = this->_kmerSize;
	Type un = 1;
	size_t _shift_val = Type::getSize() -8;
	Type kmerMask = (un << (ks*2)) - un;
	size_t shift = 2*(ks-1);
	
	/** We iterate the superkmers from the table. */
	for (it->first(); !it->isDone(); it->next())
	{
	    /** A superkmer is encoded with two successive Type objects, so we read both of them. */
		Type superk = it->item();
		it->next();
		Type seedk = it->item();
		
		compactedK =  superk;
		nbK = (compactedK >> _shift_val).getVal() & 255; // 8 bits poids fort = cpt //todo for large k values
		rem = nbK;
		
		Type temp = seedk;
		Type rev_temp = revcomp(temp,ks);
		Type newnt ;
		Type mink;
		
		/** We loop over each kmer of the current superkmer. */
		for (int ii=0; ii<nbK; ii++,rem--)
		{
			mink = std::min (rev_temp, temp);
			
			/** We insert the kmer into the hash. */
			hash.increment (mink);
							
			if(rem < 2) break;
			newnt =  ( superk >> ( 2*(rem-2)) ) & 3 ;
			
			temp = ((temp << 2 ) |  newnt   ) & kmerMask;
			newnt =  Type(comp_NT[newnt.getVal()]) ;
			rev_temp = ((rev_temp >> 2 ) |  (newnt << shift) ) & kmerMask;
		}
	}

	/** We loop over the solid kmers map.
	 * NOTE !!! we want the items to be sorted by kmer values (see finalize part of debloom). */
	Iterator < Abundance<Type> >* itKmerAbundance = hash.iterator(true);
	LOCAL (itKmerAbundance);

	for (itKmerAbundance->first(); !itKmerAbundance->isDone(); itKmerAbundance->next())
	{
		/** We may add this kmer to the solid kmers bag. */
	   this->insert ((Count&) itKmerAbundance->item());
	}
	
	this->_progress->inc (this->_pInfo.getNbKmer(this->_parti_num) ); // this->_pInfo->getNbKmer(this->_parti_num)  kmers.size()
};

/*********************************************************************
        #     #  #######   #####   #######  #######  ######
        #     #  #        #     #     #     #     #  #     #
        #     #  #        #           #     #     #  #     #
        #     #  #####    #           #     #     #  ######
         #   #   #        #           #     #     #  #   #
          # #    #        #     #     #     #     #  #    #
           #     #######   #####      #     #######  #     #
*********************************************************************/

template<size_t span>
class SuperKReader
{
	typedef typename Kmer<span>::Type  Type;
public:

	void operator() (Type& elem)
	{
		//reading elems by pairs
		if(_first)
		{
			_superk = elem;
			_first = false;
		}
		else
		{
			_seedk = elem;
			
			Type compactedK;
			
			compactedK =  _superk;
			u_int8_t nbK = (compactedK >> _shift_val).getVal()  & 255; // 8 bits poids fort = cpt
			u_int8_t rem = nbK;
			
			Type temp = _seedk;
			Type rev_temp = revcomp(temp,_kmerSize);
			Type newnt ;
			Type mink, prev_mink;
			uint64_t idx;
			
			bool prev_which =  (temp < rev_temp );
			int kx_size = -1; //next loop start at ii=0, first kmer will put it at 0
			Type radix_kxmer_forward =  (temp & _mask_radix) >> ((_kmerSize - 4)*2);
			Type  first_revk, kinsert,radix_kxmer;
			
			if(!prev_which) first_revk = rev_temp;
			
			u_int8_t rid;

			for (int ii=0; ii< nbK; ii++,rem--)
			{
				bool which =  (temp < rev_temp );
				mink = which ? temp : rev_temp;
				
				if (which != prev_which || kx_size >= _kx) // kxmer_size = 1
				{
					//output kxmer size kx_size,radix_kxmer
					//kx mer is composed of superKp[ii-1] superKp[ii-2] .. superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
					
					if(prev_which)
					{
						radix_kxmer = radix_kxmer_forward;
						kinsert = prev_mink;
					}
					else // si revcomp, le radix du kxmer est le debut du dernier kmer
					{
						//previous mink
						radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
						kinsert = first_revk;
					}
					
					//record kxmer
					rid = radix_kxmer.getVal();
					//idx = _r_idx[IX(kx_size,rid)]++;
					idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
					
					_radix_kmers [IX(kx_size,rid)][ idx] = kinsert << ((4-kx_size)*2);  //[kx_size][rid]
					if (_bankIdMatrix)  { _bankIdMatrix[IX(kx_size,rid)][ idx] = _bankId; }
					
					radix_kxmer_forward =  (mink & _mask_radix) >> _shift_radix;
					kx_size =0;
					
					if(!which) first_revk = rev_temp;
				}
				else
				{
					kx_size++;
				}
			
				prev_which = which ;
				prev_mink = mink;
				
				if(rem < 2) break; //no more kmers in this superkmer, the last one has just been eaten
				newnt =  ( _superk >> ( 2*(rem-2)) ) & 3 ;
				
				temp = ((temp << 2 ) |  newnt   ) & _kmerMask;
				newnt =  Type(comp_NT[newnt.getVal()]) ;
				rev_temp = ((rev_temp >> 2 ) |  (newnt << _shift) ) & _kmerMask;
			}
			
			//record last kxmer prev_mink et monk ?
			if(prev_which)
			{
				radix_kxmer = radix_kxmer_forward;
				kinsert = prev_mink;
			}
			else // si revcomp, le radix du kxmer est le debut du dernier kmer
			{
				//previous mink
				radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
				kinsert = first_revk;
			}
			
			//record kxmer
			rid = radix_kxmer.getVal();
			//idx = _r_idx[IX(kx_size,rid)]++;
			idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
					
			_radix_kmers [IX(kx_size,rid)][ idx] = kinsert << ((4-kx_size)*2);   // [kx_size][rid]
			if (_bankIdMatrix)  { _bankIdMatrix[IX(kx_size,rid)][ idx] = _bankId; }

			_first = true;
		}
	}
	
	SuperKReader (size_t kmerSize,  uint64_t * r_idx, Type** radix_kmers, u_int8_t** bankIdMatrix, size_t bankId=0)
	: _first(true) ,_kmerSize (kmerSize), _r_idx (r_idx), _radix_kmers(radix_kmers), _bankIdMatrix(bankIdMatrix), _kx(4), _bankId(bankId)
	 {
		 Type un = 1;
		 _kmerMask    = (un << (kmerSize*2)) - un;
		 _mask_radix  = Type((int64_t) 255);
		 _mask_radix  = _mask_radix << ((_kmerSize - 4)*2);
		 _shift       = 2*(kmerSize-1);
		 _shift_val   = un.getSize() -8;
		 _shift_radix = ((kmerSize - 4)*2); // radix is 4 nt long
	}
	
private :

	size_t _kmerSize;
	size_t _shift ;
	size_t _shift_val ;
	size_t _shift_radix ;
	int    _kx;
	Type** _radix_kmers;
	u_int8_t** _bankIdMatrix;
	uint64_t* _r_idx ;
	bool _first;
	Type _superk, _seedk;
	Type _radix, _mask_radix ;
	Type _kmerMask;
	size_t _bankId;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
/** in this scheme we count k-mers in a partition by sorting a vector*/
template<size_t span>
PartitionsByVectorCommand<span>:: PartitionsByVectorCommand (
    Bag<Count>*         solidKmers,
    Iterable<Type>&     partition,
    IHistogram*         histogram,
    ISynchronizer*      synchro,
    u_int64_t&          totalKmerNbRef,
    pair<size_t,size_t> abundance,
    IteratorListener*   progress,
    TimeInfo&           timeInfo,
    PartiInfo<5>&       pInfo,
    int                 parti,
    size_t              nbCores,
    size_t              kmerSize,
    MemAllocator&       pool,
    size_t              cacheSize,
    KmerSolidityKind    solidityKind,
    vector<size_t>&     offsets
)
    : PartitionsCommand<span> (
        solidKmers, partition, histogram, synchro, totalKmerNbRef, abundance, progress, timeInfo, pInfo,parti,nbCores,kmerSize,pool,cacheSize),
        _radix_kmers (0), _bankIdMatrix(0), _radix_sizes(0), _r_idx(0), _solidityKind(solidityKind), _nbItemsPerBankPerPart(offsets)
{
    _dispatcher = new Dispatcher (this->_nbCores);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
PartitionsByVectorCommand<span>:: ~PartitionsByVectorCommand ()
{
    if (_dispatcher)  { delete _dispatcher; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsByVectorCommand<span>::execute ()
{
    /** We check that we got something. */
    if (this->_partition.getNbItems() == 0)  {  return;  }

    /** We configure tables. */
    _radix_kmers  = (Type**)     MALLOC (256*(KX+1)*sizeof(Type*)); //make the first dims static ?  5*256
    _radix_sizes  = (uint64_t*)  MALLOC (256*(KX+1)*sizeof(uint64_t));
    _r_idx        = (uint64_t*)  CALLOC (256*(KX+1),sizeof(uint64_t));

    /** We need extra information for kmers counting in case of several input banks. */
    if (_nbItemsPerBankPerPart.size() > 1) { _bankIdMatrix = (u_int8_t**) MALLOC (256*(KX+1)*sizeof(u_int8_t*)); }
    else                                   { _bankIdMatrix = 0; }

    /** We have 3 phases here: read, sort and dump. */
    executeRead ();
    executeSort ();
    executeDump ();

    /** We cleanup tables. */
    FREE (_radix_sizes) ;
    FREE (_radix_kmers);
    FREE (_r_idx);
    if (_bankIdMatrix)  { FREE (_bankIdMatrix); }

    /** We update the progress bar. */
    this->_progress->inc (this->_pInfo.getNbKmer(this->_parti_num) );
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsByVectorCommand<span>::executeRead ()
{
    TIME_INFO (this->_timeInfo, "1.read");

    /** Recall that the attribute _offsets has a size equals to the number of banks + 1 as input
     * and for each bank, it holds the number of items found for the currently processed partition.
     *
     *               bank0   bank1   ...   bankI
     *   offsets :    xxx     xxx           xxx
     *               <------------------------->
     *                current partition content
     */
    DEBUG (("_offsets.size=%d  OFFSETS: ", _nbItemsPerBankPerPart.size() ));
    for (size_t j=0; j<_nbItemsPerBankPerPart.size(); j++)  {  DEBUG (("%6d ", _nbItemsPerBankPerPart[j]));  }  DEBUG (("\n"));

    uint64_t sum_nbxmer =0;
    for (int xx=0; xx< (KX+1); xx++)
    {
        for (int ii=0; ii< 256; ii++)
        {
            /** Shortcut. */
            size_t nbKmers = this->_pInfo.getNbKmer(this->_parti_num,ii,xx);

            //use memory pool here to avoid memory fragmentation
            _radix_kmers  [IX(xx,ii)] = (Type*)     this->_pool.pool_malloc (nbKmers * sizeof(Type),     "kmers alloc");
            _radix_sizes  [IX(xx,ii)] = nbKmers;

            if (_bankIdMatrix)
            {
                _bankIdMatrix [IX(xx,ii)] = (u_int8_t*) this->_pool.pool_malloc (nbKmers * sizeof(u_int8_t), "bank ids alloc");
            }

            sum_nbxmer +=  nbKmers;
        }
    }

    DEBUG (("PartitionsByVectorCommand<span>::executeRead:  fillsolid parti num %i  by vector  nb kxmer / nbkmers      %lli / %lli     %f   with %zu nbcores \n",
        this->_parti_num, sum_nbxmer, this->_pInfo.getNbKmer(this->_parti_num),
        (double) sum_nbxmer /  this->_pInfo.getNbKmer(this->_parti_num),this->_nbCores
    ));

    /** HOW TO COUNT KMERS BY SET OF READS ?
     * Now, we are going to read the temporary partition built during the previous phase and fill
     * the _radix_kmers attribute. We also need to know in _radix_kmers what is the contribution of
     * each bank. We therefore need to iterate the current partition by bank (using the information
     * of _offsets). */

    if (_bankIdMatrix)
    {
        /** We create an iterator over all the items. */
        Iterator<Type>* itGlobal = this->_partition.iterator();
        LOCAL (itGlobal);

        /** We iterate the banks. */
        for (size_t b=0; b<_nbItemsPerBankPerPart.size(); b++)
        {
            /** We truncate the global iterator.
             * NB : we initialize (ie call 'first') the global iterator only at first call (for b==0). */
            Iterator<Type>* itLocal = new TruncateIterator<Type> (*itGlobal, _nbItemsPerBankPerPart[b], b==0 ? true : false);
            LOCAL (itLocal);

            /** We iterate this local iterator. */
            _dispatcher->iterate (itLocal, SuperKReader<span>  (this->_kmerSize, _r_idx, _radix_kmers, _bankIdMatrix, b), 10000); //must be even , reading by pairs
        }

        /** We check that the global iterator is finished. */
        if (itGlobal->isDone() == false)  { throw Exception ("PartitionsByVectorCommand: iteration should be finished"); }
    }
    else
    {
        /** We iterate the superkmers. */
        _dispatcher->iterate (this->_partition.iterator(), SuperKReader<span>  (this->_kmerSize, _r_idx, _radix_kmers, 0, 0), 10000); //must be even , reading by pairs
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
class SortCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
public:
    typedef typename Kmer<span>::Type  Type;

    /** Constructor. */
    SortCommand (Type** kmervec, u_int8_t** bankIdMatrix, int begin, int end, uint64_t* radix_sizes)
        : _radix_kmers(kmervec), _bankIdMatrix(bankIdMatrix), _deb(begin), _fin(end), _radix_sizes(radix_sizes) {}

    /** */
    void execute ()
    {
        vector<size_t> idx;
        vector<Tmp>    tmp;

        for (int ii=_deb; ii <=_fin; ii++)
        {
            if (_radix_sizes[ii] > 0)
            {
                /** Shortcuts. */
                Type* kmers = _radix_kmers  [ii];

                if (_bankIdMatrix)
                {
                    /** NOT OPTIMAL AT ALL... in particular we have to use 'idx' and 'tmp' vectors
                     * which may use (a lot of ?) memory. */

                    /** Shortcut. */
                    u_int8_t* banksId = _bankIdMatrix [ii];

                    /** NOTE: we sort the indexes, not the items. */
                    idx.resize (_radix_sizes[ii]);
                    for (int i=0; i<idx.size(); i++)  { idx[i]=i; }

                    std::sort (idx.begin(), idx.end(), Cmp(kmers));

                    /** Now, we have to reorder the two provided vectors with the same order. */
                    tmp.resize (idx.size());
                    for (int i=0; i<idx.size(); i++)
                    {
                        tmp[i].kmer = kmers  [idx[i]];
                        tmp[i].id   = banksId[idx[i]];
                    }
                    for (int i=0; i<idx.size(); i++)
                    {
                        kmers  [i] = tmp[i].kmer;
                        banksId[i] = tmp[i].id;
                    }
                }
                else
                {
                    std::sort (&kmers[0] , &kmers[ _radix_sizes[ii]]);
                }
            }
        }
    }

private :

    struct Tmp { Type kmer;  u_int8_t id;};

    struct Cmp
    {
        Type* _kmers;
        Cmp (Type* kmers) : _kmers(kmers) {}
        bool operator() (size_t a, size_t b)  { return _kmers[a] < _kmers[b]; }
    };

    int        _deb;
    int        _fin;
    Type**     _radix_kmers;
    u_int8_t** _bankIdMatrix;
    uint64_t*  _radix_sizes;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsByVectorCommand<span>::executeSort ()
{
    TIME_INFO (this->_timeInfo, "2.sort");

    vector<ICommand*> cmds;

    int nwork = 256 / this->_nbCores;

    for (int xx=0; xx < (KX+1); xx++)
    {
        cmds.clear();

        //fill cmd work vector
        for (int tid=0; tid < this->_nbCores; tid++)
        {
            int deb = 0 + tid * nwork;
            int fin = (tid+1) * nwork -1; // thread will do inclusive range [begin -- end ]
            if(tid == this->_nbCores-1)  { fin = 255; }

            // mettre dans le  SortCommand le master radix_kmers et range a traiter
            cmds.push_back (new SortCommand<span> (
                _radix_kmers+ IX(xx,0),
                (_bankIdMatrix ? _bankIdMatrix+ IX(xx,0) : 0),
                deb, fin,
                _radix_sizes + IX(xx,0)
            ));
        }

        _dispatcher->dispatchCommands (cmds, 0);
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
class KxmerPointer
{
public:
    typedef typename Kmer<span>::Type  Type;

    //on lui passe le vector dun kxmer //std::vector<vector<Type> >  & kmervec
    KxmerPointer (
        Type**      kmervec,
        int         prefix_size,
        int         x_size,
        int         min_radix,
        int         max_radix,
        int         kmerSize,
        uint64_t*   radix_sizes,
        u_int8_t**  bankIdMatrix
    )
        : _kxmers(kmervec), _prefix_size(prefix_size), _x_size(x_size),
        _cur_idx(-1) ,_kmerSize(kmerSize),_low_radix(min_radix),_high_radix(max_radix), _radix_sizes(radix_sizes), _bankIdMatrix(0)
    {
        _idx_radix = min_radix;
        Type un = 1;
        _kmerMask = (un << (_kmerSize*2)) - un;

        _shift_size = ( (4 - _prefix_size) *2) ;
        _radixMask = Type(_idx_radix) ;
        _radixMask = _radixMask << ((_kmerSize-4)*2);
        _radixMask = _radixMask  << (2*_prefix_size)  ;

        if (bankIdMatrix) { _bankIdMatrix = bankIdMatrix + IX(x_size,0); }
    }

    /** */
    inline bool next ()
    {
        _cur_idx++;

        // go to next non empty radix
        while(_idx_radix<= _high_radix &&_cur_idx >=   _radix_sizes[_idx_radix])
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

    /** */
    inline Type    value   () const  {  return ( ((_kxmers[_idx_radix][_cur_idx]) >> _shift_size)  |  _radixMask  ) & _kmerMask ;  }

    /** */
    inline u_int8_t getBankId () const  {  return _bankIdMatrix ? _bankIdMatrix [_idx_radix][_cur_idx] : 0;  }

private :

    Type**      _kxmers;
    u_int8_t**  _bankIdMatrix;
    uint64_t*   _radix_sizes;
    int64_t     _cur_idx;
    Type        _cur_radix;
    Type        _kmerMask;
    Type        _radixMask;
    int         _idx_radix;
    int         _low_radix;
    int         _high_radix;
    int         _shift_size;
    int         _prefix_size;
    int         _kmerSize;
    int         _x_size; //x size of the _kxmersarray
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/

template<size_t span>
void PartitionsByVectorCommand<span>::executeDump ()
{
    TIME_INFO (this->_timeInfo, "3.dump");

    int nbkxpointers = 453; //6 for k1 mer, 27 for k2mer, 112 for k3mer  453 for k4mer
    vector< KxmerPointer<span>*> vec_pointer (nbkxpointers);
    int best_p;

    std::priority_queue< kxp, std::vector<kxp>,kxpcomp > pq;

    SolidityCounter solidCounter (_solidityKind, this->_abundance, _nbItemsPerBankPerPart.size());

    Type previous_kmer ;

    //init the pointers to the 6 arrays
    int pidx =0;

    ////////////////////////////////////////////////
    ////-------------k0 pointers-----------/////////
    ////////////////////////////////////////////////

    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(0,0) ,0,0,0,255,this->_kmerSize, _radix_sizes + IX(0,0), _bankIdMatrix); // vec, prefix size, kxsize , radix min, radix max ,ksize

    ////////////////////////////////////////////////
    ////-------------k1 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(1,0) ,0,1,0,255,this->_kmerSize, _radix_sizes + IX(1, 0), _bankIdMatrix);
    int lowr = 0;
    int maxr = 63;

    //prefix1
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(1,0) ,1,1,lowr,maxr,this->_kmerSize, _radix_sizes + IX(1, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    ////////////////////////////////////////////////
    ////-------------k2 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),0,2,0,255,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);

    //prefix1
    lowr = 0; maxr = 63;
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),1,2,lowr,maxr,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    //prefix2
    lowr = 0; maxr = 15;
    for(unsigned int ii=0; ii<16; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),2,2,lowr,maxr,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
        lowr += 16;
        maxr += 16;
    }

    ////////////////////////////////////////////////
    ////-------------k3 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),0,3,0,255,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);

    //prefix1
    lowr = 0; maxr = 63;
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),1,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    //prefix2
    lowr = 0; maxr = 15;
    for(unsigned int ii=0; ii<16; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),2,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
        lowr += 16;
        maxr += 16;
    }

    //prefix3
    lowr = 0; maxr = 3;
    for(unsigned int ii=0; ii<64; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),3,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
        lowr += 4;
        maxr += 4;
    }

    ////////////////////////////////////////////////
    ////-------------k4 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),0,4,0,255,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);

    //prefix1
    lowr = 0; maxr = 63;
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),1,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    //prefix2
    lowr = 0; maxr = 15;
    for(unsigned int ii=0; ii<16; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),2,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 16;
        maxr += 16;
    }

    //prefix3
    lowr = 0; maxr = 3;
    for(unsigned int ii=0; ii<64; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),3,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 4;
        maxr += 4;
    }

    //prefix4
    lowr = 0; maxr = 0;
    for(unsigned int ii=0; ii<256; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),4,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 1;
        maxr += 1;
    }

    //fill the  priority queue with the first elems
    for (int ii=0; ii<nbkxpointers; ii++)
    {
        if(vec_pointer[ii]->next())  {  pq.push(kxp(ii,vec_pointer[ii]->value()));  }
    }

    if (pq.size() != 0) // everything empty, no kmer at all
    {
        //get first pointer
        best_p = pq.top().first ; pq.pop();

        previous_kmer = vec_pointer[best_p]->value();

        solidCounter.init (vec_pointer[best_p]->getBankId());

        //merge-scan all 'virtual' arrays and output counts
        while (1)
        {
            //go forward in this array or in new array of reaches end of this one
            if (! vec_pointer[best_p]->next())
            {
                //reaches end of one array
                if(pq.size() == 0) break; //everything done

                //otherwise get new best
                best_p = pq.top().first ; pq.pop();
            }

            if (vec_pointer[best_p]->value() != previous_kmer )
            {
                //if diff, changes to new array, get new min pointer
                pq.push(kxp(best_p,vec_pointer[best_p]->value())); //push new val of this pointer in pq, will be counted later

                best_p = pq.top().first ; pq.pop();

                //if new best is diff, this is the end of this kmer
                if(vec_pointer[best_p]->value()!=previous_kmer )
                {
                    this->insert (previous_kmer, solidCounter);

                    solidCounter.init (vec_pointer[best_p]->getBankId());
                    previous_kmer = vec_pointer[best_p]->value();
                }
                else
                {
                    solidCounter.increase (vec_pointer[best_p]->getBankId());
                }
            }
            else
            {
                solidCounter.increase (vec_pointer[best_p]->getBankId());
            }
        }

        //last elem
        this->insert (previous_kmer, solidCounter);
    }

    /** Cleanup. */
    for (size_t ii=0; ii<nbkxpointers; ii++)  {  delete vec_pointer[ii];  }
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class PartitionsCommand <KSIZE_1>;
template class PartitionsCommand <KSIZE_2>;
template class PartitionsCommand <KSIZE_3>;
template class PartitionsCommand <KSIZE_4>;

template class PartitionsByHashCommand <KSIZE_1>;
template class PartitionsByHashCommand <KSIZE_2>;
template class PartitionsByHashCommand <KSIZE_3>;
template class PartitionsByHashCommand <KSIZE_4>;

template class PartitionsByVectorCommand <KSIZE_1>;
template class PartitionsByVectorCommand <KSIZE_2>;
template class PartitionsByVectorCommand <KSIZE_3>;
template class PartitionsByVectorCommand <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
