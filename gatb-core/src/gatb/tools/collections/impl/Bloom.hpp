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

/** \file Bloom.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bloom implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/system/api/types.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <bitset>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

extern u_int8_t bit_mask [];

/********************************************************************************/
/**
 */
template <typename Item> class HashFunctors
{
public:

    /** */
    HashFunctors (size_t nbFct, u_int64_t seed=0) : _nbFct(nbFct), user_seed(seed)
    {
        generate_hash_seed ();
    }

    /** */
    u_int64_t operator ()  (const Item& key, size_t idx)  {  return hash1 (key, seed_tab[idx]);  }

private:

    /** */
    void generate_hash_seed ()
    {
        static const u_int64_t rbase[NSEEDSBLOOM] =
        {
            0xAAAAAAAA55555555ULL,  0x33333333CCCCCCCCULL,  0x6666666699999999ULL,  0xB5B5B5B54B4B4B4BULL,
            0xAA55AA5555335533ULL,  0x33CC33CCCC66CC66ULL,  0x6699669999B599B5ULL,  0xB54BB54B4BAA4BAAULL,
            0xAA33AA3355CC55CCULL,  0x33663366CC99CC99ULL
        };

        for (size_t i=0; i<NSEEDSBLOOM; ++i)  {  seed_tab[i] = rbase[i];  }
        for (size_t i=0; i<NSEEDSBLOOM; ++i)  {  seed_tab[i] = seed_tab[i] * seed_tab[(i+3) % NSEEDSBLOOM] + user_seed ;  }
    }

    size_t _nbFct;

    static const size_t NSEEDSBLOOM = 10;
    u_int64_t seed_tab[NSEEDSBLOOM];
    u_int64_t user_seed;
};

/********************************************************************************/

/** \brief Bloom interface
 *
 * We define a Bloom filter has a container (something that tells whether an item is here or not) and
 * a bag (something we can put items into).
 *
 * As expected, there is no Iterable interface here.
 */
template <typename Item> class IBloom : public Container<Item>, public Bag<Item>, public system::SmartPointer
{
public:
    virtual ~IBloom() {}

    virtual u_int8_t*& getArray    () = 0;
    virtual u_int64_t  getSize     () = 0;
    virtual u_int64_t  getBitSize  () = 0;
    virtual size_t     getNbHash   () const = 0;

	virtual std::bitset<4> contains4 (const Item& item, bool right) = 0;

    virtual std::bitset<8> contains8 (const Item& item) = 0;

    virtual std::string  getName   () const  = 0;
    virtual unsigned long  weight () = 0;
};

/********************************************************************************/

/** \brief Bloom filter implementation
 */
template <typename Item> class BloomContainer : public IBloom<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use */
    BloomContainer (u_int64_t tai_bloom, size_t nbHash = 4)
        : _hash(nbHash), n_hash_func(nbHash), blooma(0), tai(tai_bloom), nchar(0), isSizePowOf2(false)
    {
        nchar  = (1+tai/8LL);
        blooma = (unsigned char *) MALLOC (nchar*sizeof(unsigned char)); // 1 bit per elem
        system::impl::System::memory().memset (blooma, 0, nchar*sizeof(unsigned char));

        /** We look whether the provided size is a power of 2 or not.
         *   => if we have a power of two, we can optimize the modulo operations. */
        isSizePowOf2 = (tai && !(tai & (tai - 1)));

        /** In case we have a power of 2^N, we set the size to 2^N-1 and use the classical trick:
         *     a % 2^N  <=>  a & (2^N-1)
         * */
        if (isSizePowOf2)  {  tai --;  }
    }

    /** Destructor. */
    virtual ~BloomContainer ()
    {
        system::impl::System::memory().free (blooma);
    }

    /** */
    size_t getNbHash () const { return n_hash_func; }

    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        if (isSizePowOf2)
        {
            for (size_t i=0; i<n_hash_func; i++)
            {
                u_int64_t h1 = _hash (item,i) & tai;
               // if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
				if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }

            }
        }
        else
        {
            for (size_t i=0; i<n_hash_func; i++)
            {
                u_int64_t h1 = _hash (item,i) % tai;
               // if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
				if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }

            }
        }
        return true;
    }

	virtual std::bitset<4> contains4 (const Item& item, bool right)  { return std::bitset<4>();} //not implem by default

    virtual std::bitset<8> contains8 (const Item& item)  { return std::bitset<8>(); }

    /** */
    virtual u_int8_t*& getArray     ()  { return blooma; }
    virtual u_int64_t  getSize      ()  { return nchar;  }
    virtual u_int64_t  getBitSize   ()  { return tai;    }

    virtual std::string  getName    () const  = 0;

protected:

    HashFunctors<Item> _hash;
    int n_hash_func;

    u_int8_t* blooma;
    u_int64_t tai;
    u_int64_t nchar;
    bool      isSizePowOf2;
};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class Bloom : public BloomContainer<Item>
{
public:

    /** \copydoc BloomContainer::BloomContainer */
    Bloom (u_int64_t tai_bloom, size_t nbHash = 4)  : BloomContainer<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        if (this->isSizePowOf2)
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) & this->tai;
                this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
            }
        }
        else
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) % this->tai;
                this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
            }
        }
    }

    /** \copydoc Bag::flush */
    void flush ()  {}

    /** */
    std::string  getName () const { return "Bloom"; }

    /** */
    void dump (const char* filename)
    {
        FILE* file = fopen(filename,"wb");
        fwrite (this->blooma, sizeof(unsigned char), this->nchar, file);
        fclose (file);
    }

    // return the number of 1's in the Bloom, nibble by nibble
    unsigned long weight()
    {
        const unsigned char oneBits[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
        long weight = 0;
        for(uint64_t index = 0; index < this->nchar; index++)
        {
            unsigned char current_char = this->blooma[index];
            weight += oneBits[current_char&0x0f];
            weight += oneBits[current_char>>4];
        }
        return weight;
    }



};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class BloomNull : public IBloom<Item>
{
public:

    virtual ~BloomNull() {}

    u_int8_t*& getArray    () { return a; }
    u_int64_t  getSize     () { return 0; }
    u_int64_t  getBitSize  () { return 0; }
    size_t     getNbHash   () const { return 0; }

    virtual std::string  getName   () const  { return "BloomNull"; }

    std::bitset<4> contains4 (const Item& item, bool right)  {return std::bitset<4>();}

	std::bitset<8> contains8 (const Item& item)  { return std::bitset<8>(); }

	
    bool contains (const Item& item) { return false; }
    void insert (const Item& item) {}
    void flush ()  {}
    unsigned long weight ()  { return 0;}
private:
    u_int8_t* a;
};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class BloomSynchronized : public Bloom<Item>
{
public:

    /** \copydoc BloomContainer::BloomContainer */
    BloomSynchronized (u_int64_t tai_bloom, size_t nbHash = 4)  : Bloom<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        if (this->isSizePowOf2)
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) & this->tai;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        }
        else
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) % this->tai;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        }
    }

    /** */
    std::string  getName () const { return "basic"; }
};
    

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class BloomCacheCoherent : public Bloom<Item>
{
public:
    
    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use
     * \param[in] block_nbits : size of the block (actual 2^nbits) */
    BloomCacheCoherent (u_int64_t tai_bloom, size_t nbHash = 4,size_t block_nbits = 12)  : Bloom<Item> (tai_bloom + 2*(1<<block_nbits), nbHash),_nbits_BlockSize(block_nbits)
    {
        _mask_block = (1<<_nbits_BlockSize) - 1;
        _reduced_tai = this->tai -  2*(1<<_nbits_BlockSize) ;//2* for neighbor coherent
    }
    
    
    //for insert, no prefetch, perf is not important
     /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {

        u_int64_t h0;

        h0 = this->_hash (item,0) % _reduced_tai;
        
        __sync_fetch_and_or (this->blooma + (h0 >> 3), bit_mask[h0 & 7]);

            for (size_t i=1; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = h0  + (simplehash16( item, i) & _mask_block )   ;
			//	u_int64_t h1 = h0  +  (this->_hash (item,i)  & _mask_block );

                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        
    }
    
    /** */
    std::string  getName () const { return "cache"; }

    /** */
    u_int64_t  getBitSize   ()  { return _reduced_tai;    }
        
    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {

        u_int64_t tab_keys [20];
        u_int64_t h0;

        h0 = this->_hash (item,0) % _reduced_tai;
        __builtin_prefetch(&(this->blooma [h0 >> 3] ), 0, 3); //preparing for read

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
           tab_keys[i] =  h0  + (simplehash16( item, i) & _mask_block );// with simplest hash
		//	tab_keys[i] =  h0  + (this->_hash (item,i) & _mask_block );// with simplest hash

        }
        
        
//        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) != bit_mask[h0 & 7])  {  return false;  }
        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) ==0 )  {  return false;  }
        
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = tab_keys[i];
           // if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
			if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }

        }
        return true;
        
    }
    
    unsigned long weight()
    {
        std::cout << "GATB Warning, called unimpplemented BloomCacheCoherent::weight()\n";
        return 0; // not implemented
    }
    
protected:
    u_int64_t _mask_block;
    size_t _nbits_BlockSize;
    u_int64_t _reduced_tai;
};

	
/********************************************************************************/
	
template <typename Item> class BloomNeighborCoherent : public BloomCacheCoherent<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use
     * \param[in] block_nbits : size of the block (actual 2^nbits) */
    BloomNeighborCoherent (u_int64_t tai_bloom, size_t kmersize , size_t nbHash = 4,size_t block_nbits = 12 )  :
    BloomCacheCoherent<Item> (tai_bloom , nbHash,block_nbits), _kmerSize(kmersize)
    {
        cano2[ 0] = 0;
        cano2[ 1] = 1;
        cano2[ 2] = 2;
        cano2[ 3] = 3;
        cano2[ 4] = 4;
        cano2[ 5] = 5;
        cano2[ 6] = 3;
        cano2[ 7] = 7;
        cano2[ 8] = 8;
        cano2[ 9] = 9;
        cano2[10] = 0;
        cano2[11] = 4;
        cano2[12] = 9;
        cano2[13] = 13;
        cano2[14] = 1;
        cano2[15] = 5;

        Item un = 1;
        _maskkm2 = (un << ((_kmerSize-2)*2)) - un;
        _kmerMask = (un << (_kmerSize*2)) - un;

        Item trois = 3;

        _prefmask = trois << ((_kmerSize-1)*2); //bug was here 3 instead of item trois
    }

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        u_int64_t h0;
        u_int64_t racine;

        Item suffix = item & 3 ;
        Item prefix = (item & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        prefix = prefix  & 15 ;

        u_int64_t pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        Item hashpart = ( item >> 2 ) & _maskkm2 ;  // delete 1 nt at each side
        Item rev =  revcomp(hashpart,_kmerSize-2);
        if(rev<hashpart) hashpart = rev; //transform to canonical

//			Item km = item;
//			rev =  revcomp(km,_kmerSize);
//			if(rev < km) km = rev; //transform to canonical
        
        racine = ((this->_hash (hashpart,0) ) % this->_reduced_tai) ;
        //h0 = ((this->_hash (item >> 2,0) ) % this->_reduced_tai)  + (suffix_val & this->_mask_block);
        //h0 = racine + (this->_hash (km,0)  & this->_mask_block);
        h0 = racine + (pref_val );
        __sync_fetch_and_or (this->blooma + (h0 >> 3), bit_mask[h0 & 7]);

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = h0  + ( (simplehash16( hashpart, i))  & this->_mask_block )   ;
            //	u_int64_t h1 = racine  + ( (simplehash16( km, i))  & this->_mask_block )   ; //ceci avec simplehash+8  semble ok
            //	u_int64_t h1 = h0  +  ( (this->_hash (item>>2,i)+ suffix_val)  & _mask_block );
            __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
        }
    }

    /** */
    std::string  getName () const { return "neighbor"; }

    /** */
    u_int64_t  getBitSize   ()  { return this->_reduced_tai;    }

    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        u_int64_t racine;

        Item suffix = item & 3 ;
        Item prefix = (item & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        prefix = prefix  & 15 ;

        u_int64_t pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        Item hashpart = ( item >> 2 ) & _maskkm2 ;  // delete 1 nt at each side
        Item rev =  revcomp(hashpart,_kmerSize-2);
        if(rev<hashpart) hashpart = rev; //transform to canonical

        // Item km = item;
        // rev =  revcomp(km,_kmerSize);
        // if(rev < km) km = rev; //transform to canonical

        u_int64_t tab_keys [20];
        u_int64_t h0;

        racine = ((this->_hash (hashpart,0) ) % this->_reduced_tai) ;
        //h0 = racine + (this->_hash (km,0)  & this->_mask_block);
        h0 = racine + (pref_val  );
        //printf("h0 %llu\n",h0);

        __builtin_prefetch(&(this->blooma [h0 >> 3] ), 0, 3); //preparing for read

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            tab_keys[i] =  h0  + (  (simplehash16( hashpart, i)  ) & this->_mask_block );// with simplest hash
            // tab_keys[i] =  racine  + (  (simplehash16( km, i)  ) & this->_mask_block );
            // tab_keys[i] =  h0  + (  (this->_hash (item>>2,i) + suffix_val ) & _mask_block );// with simplest hash
        }

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0 )  {  return false;  } //was != bit_mask[h0 & 7]

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = tab_keys[i];
            if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0 )  {  return false;  } //was != bit_mask[h0 & 7]
        }
        return true;
    }

    //ask for all 4 neighbors of item in one call
    // give it CAAA
    // if right == true  it wil test the 4 neighbors AAAA, AAAC, AAAT, AAAG
    // if right == false : left extension   ACAA, CCAA , TCAA  , GCAA

    std::bitset<4> contains4 (const Item& item, bool right)
    {
        u_int64_t h0, i0, j0, k0;
        Item elem,hashpart,rev ;
        Item un = 1;
        Item deux = 2;
        Item trois = 3;

        size_t shifts = (_kmerSize -1)*2;

        if (right)  {  elem = (item << 2) & _kmerMask ;  }
        else        {  elem = (item >> 2) ;              }

        //get the canonical of middle part
        hashpart = ( elem >> 2 ) & _maskkm2 ;
        rev =  revcomp(hashpart,_kmerSize-2);
        if(rev<hashpart) hashpart = rev;

        u_int64_t racine = ((this->_hash (hashpart,0) ) % this->_reduced_tai) ;

        __builtin_prefetch(&(this->blooma [racine >> 3] ), 0, 3); //preparing for read

        Item tmp,suffix,prefix;
        u_int64_t pref_val;

        //with val of prefix+suffix  for neighbor shift

        tmp = elem;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        h0 = racine + (pref_val  & this->_mask_block);

        if(right) tmp = elem+un;
        else tmp = elem + (un<<shifts) ;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        i0 = racine + (pref_val  & this->_mask_block);

        if(right) tmp = elem+deux;
        else tmp = elem + (deux<<shifts) ;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        j0 = racine + (pref_val  & this->_mask_block);

        if(right) tmp = elem+trois;
        else tmp = elem + (trois<<shifts) ;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        k0 = racine + (pref_val  & this->_mask_block);

        /*
         //with full hash of kmer for neighbor shift
        tmp = elem;
        rev =  revcomp(tmp,_kmerSize);
        if(rev < tmp) tmp = rev;
        h0 = racine + (this->_hash (tmp,0)  & this->_mask_block);

        if(right) tmp = elem+un;
        else tmp = elem + (un<<shifts) ;
        rev =  revcomp(tmp,_kmerSize); // many revcomp, optim possible
        if(rev < tmp) tmp = rev;
        i0 = racine + (this->_hash (tmp,0)  & this->_mask_block);

        if(right) tmp = elem+deux;
        else tmp = elem + (deux<<shifts) ;
        rev =  revcomp(tmp,_kmerSize);
        if(rev < tmp) tmp = rev;
        j0 = racine + (this->_hash (tmp,0)  & this->_mask_block);

        if(right) tmp = elem+trois;
        else tmp = elem + (trois<<shifts) ;
        rev =  revcomp(tmp,_kmerSize);
        if(rev < tmp) tmp = rev;
        k0 = racine + (this->_hash (tmp,0)  & this->_mask_block);
         */

        u_int64_t tab_hashes [20];

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            tab_hashes[i] = simplehash16( hashpart, i) & this->_mask_block ;
        }

        std::bitset<4> resu;
        resu.set (0, true);
        resu.set (1, true);
        resu.set (2, true);
        resu.set (3, true);

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0)  {  resu.set (0, false);  }
        if ((this->blooma[i0 >> 3 ] & bit_mask[i0 & 7]) == 0)  {  resu.set (1, false);  }
        if ((this->blooma[j0 >> 3 ] & bit_mask[j0 & 7]) == 0)  {  resu.set (2, false);  }
        if ((this->blooma[k0 >> 3 ] & bit_mask[k0 & 7]) == 0)  {  resu.set (3, false);  }

        //plus rapide avec 4 boucles separees que une ci dessous avec test pour break
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  h0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (0, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  i0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (1, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  j0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (2, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  k0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (3, false); break; }
        }

        /*
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  h0 +  tab_hashes[i]   ;
            u_int64_t i1 =  i0 +  tab_hashes[i]   ;
            u_int64_t j1 =  j0 +  tab_hashes[i]   ;
            u_int64_t k1 =  k0 +  tab_hashes[i]   ;

            if (resu[0] && (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu[0]=false;  }  //test  resu[0] &&
            if (resu[1] &&(this->blooma[i1 >> 3 ] & bit_mask[i1 & 7]) == 0  )  {  resu[1]=false;  }
            if (resu[2] &&(this->blooma[j1 >> 3 ] & bit_mask[j1 & 7]) == 0  )  {  resu[2]=false;  }
            if (resu[3] &&(this->blooma[k1 >> 3 ] & bit_mask[k1 & 7]) == 0  )  {  resu[3]=false;  }

            if(resu[0]== false &&  resu[1]== false && resu[2]== false && resu[3]== false)
                break;
        }
         */

        return resu;
    }

    /** */
    std::bitset<8> contains8 (const Item& item)
    {
        std::bitset<4> resultRight = this->contains4 (item, true);
        std::bitset<4> resultLeft  = this->contains4 (item, false);
        std::bitset<8> result;
        size_t i=0;
        for (size_t j=0; j<4; j++)  { result.set (i++, resultRight[j]); }
        for (size_t j=0; j<4; j++)  { result.set (i++, resultLeft [j]); }
        return result;
    }

private:
    unsigned int cano2[16];
    Item _maskkm2;
    Item _prefmask;
    Item _kmerMask;
    size_t _kmerSize;
};
	
	
/********************************************************************************/

template <typename Item, size_t prec=1> class BloomGroupOld : public system::SmartPointer
{
public:

    typedef tools::math::LargeInt<prec> Result;

    /** */
    BloomGroupOld (u_int64_t size, size_t nbHash=4)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0)
    {
        _blooma = (Result*) MALLOC (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));
    }

    /** */
    BloomGroupOld (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);
    }

    /** */
    ~BloomGroupOld ()  {  FREE (_blooma); }

    /** */
    std::string getName () const { return "BloomGroupOld"; }

    /** */
    void insert (const Item& item, size_t idx)
    {
        static const Result ONE (1);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
#if 1
            this->_blooma[h1] |= (ONE << idx);
#else
            this->_blooma[h1].sync_fetch_and_or (ONE << idx);
#endif
        }
    }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

#if 1
for (size_t i=0; i<10; i++)
{
    for (size_t j=0; j<prec; j++)
    {
        printf ("%8x ", _blooma[i].value[j]);
    }
    printf ("\n");
}
#endif

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We allocate the array. */
            _blooma = (Result*) MALLOC (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        static const Result ONE (1);
        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            if ( (_blooma[h1] & (ONE << idx)) != (ONE << idx) )  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        static const Result ZERO (0);
        Result res = ~ZERO;

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            res &= _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;
};

/********************************************************************************/

template <typename Item, size_t prec=1> class BloomGroup : public system::SmartPointer
{
public:

    class Result
    {
    public:
        Result (u_int64_t v=0)  { memset (v); }

              u_int64_t& operator[] (size_t idx)       { return value[idx]; }
        const u_int64_t& operator[] (size_t idx) const { return value[idx]; }

        const u_int64_t* array () const { return value; }

        Result& operator&= (const Result& r)
        {
            for (size_t j=0; j<prec; j++)  { (*this)[j] &=  r[j]; }
            return *this;
        }

    private:
        u_int64_t value[prec];
        void memset (u_int64_t v)  {  system::impl::System::memory().memset (value, v, prec*sizeof(u_int64_t));  }

        friend class BloomGroup<Item,prec>;
    };

    /** */
    BloomGroup (u_int64_t size, u_int64_t maxMemory, size_t nbHash=4)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0)
    {
        printf ("BloomGroup:  size=%ld   sizeof(Result)=%d  maxMemory=%ld\n", size, sizeof(Result), maxMemory);
        if (_size*sizeof(Result) > maxMemory)
        {
            _size = maxMemory /sizeof (Result);
        }
        else
        {
            maxMemory = _size*sizeof(Result);
        }
        printf ("===> size=%ld   allocMemory=%ld\n", _size, sizeof(Result)*_size);

        _blooma = (Result*) MALLOC (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));
    }

    /** */
    BloomGroup (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);
    }

    /** */
    ~BloomGroup ()  {  FREE (_blooma); }

    /** */
    std::string getName () const { return "BloomGroup"; }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We allocate the array. */
            _blooma = (Result*) MALLOC (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** Insert an item in the 'idx'th Bloom filter
     * \param[in] item : item to be inserted into the filter
     * \param[in] idx : index of the filter */
    void insert (const Item& item, size_t idx)
    {
        u_int64_t q,mask;  euclidian(idx,q,mask);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;

#if 1
            this->_blooma[h1][q] |= mask;
#else
            __sync_fetch_and_or (this->_blooma[h1].value + q, mask);
#endif
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        u_int64_t q,mask;  euclidian(idx,q,mask);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            if ( (_blooma[h1][q] & mask) != mask )  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        Result res (~0);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            res &=  _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;

    void euclidian (size_t idx, u_int64_t& dividend, u_int64_t& mask) const
    {
        dividend = idx / (8*sizeof(u_int64_t));
        mask     = ((u_int64_t) 1) << (idx % (8*sizeof(u_int64_t)));
    }
};


/********************************************************************************/

template <typename Item, size_t prec=1> class BloomGroupCacheCoherent : public system::SmartPointer
{
public:

    typedef tools::math::LargeInt<prec> Result;

    /** */
    BloomGroupCacheCoherent (u_int64_t size, size_t nbHash=4, size_t block_nbits=12)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0), _nbits_BlockSize(block_nbits)
    {
        _size += (1<<_nbits_BlockSize);

        _blooma = (Result*) MALLOC (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

        _mask_block   = (1<<_nbits_BlockSize) - 1;
        _reduced_size = this->_size -  (1<<_nbits_BlockSize) ;
    }

    /** */
    BloomGroupCacheCoherent (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);

        _mask_block   = (1<<_nbits_BlockSize) - 1;
        _reduced_size = this->_size -  (1<<_nbits_BlockSize) ;
    }

    /** */
    ~BloomGroupCacheCoherent ()  {  system::impl::System::memory().free (_blooma); }

    /** */
    std::string getName () const { return "BloomGroupCacheCoherent"; }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write block size information. */
            file->fwrite (&_nbits_BlockSize, sizeof(_nbits_BlockSize), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We read block size information. */
            file->fread (&_nbits_BlockSize, sizeof(_nbits_BlockSize), 1);

            /** We allocate the array. */
            _blooma = (Result*) MALLOC (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void insert (const Item& item, size_t idx)
    {
        static const Result ONE (1);

        /** First hash. */
        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        this->_blooma[h0] |= (ONE << idx);

        for (size_t i=1; i<this->_nbHash; i++)
        {
            /** Other hash. */
            u_int64_t h1 = h0  + (simplehash16 (item,i) & _mask_block);
            this->_blooma[h1] |= (ONE << idx);
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        static const Result ONE  (1);

        /** First hash. */
        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        if ((this->_blooma[h0] & (ONE << idx)) != (ONE << idx))  {  return false;  }

        for (size_t i=1; i<this->_nbHash; i++)
        {
            /** Other hash. */
            u_int64_t h1 = h0  + (simplehash16 (item,i) & _mask_block);
            if ((this->_blooma[h1] & (ONE << idx)) != (ONE << idx))  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        static const Result ZERO (0);
        Result res = ~ZERO;

        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        res &= _blooma [h0];

        for (size_t i=1; i<this->_nbHash; i++)
        {
            u_int64_t h1 = h0 + (simplehash16(item,i) & _mask_block);
            res &= _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;

    u_int64_t  _mask_block;
    size_t     _nbits_BlockSize;
    u_int64_t  _reduced_size;

};

/********************************************************************************/

/** */
class BloomFactory
{
public:

    /** */
    static BloomFactory& singleton()  { static BloomFactory instance; return instance; }

    /** */
    template<typename T> IBloom<T>* createBloom (tools::misc::BloomKind kind, u_int64_t tai_bloom, size_t nbHash, size_t kmersize)
    {
        switch (kind)
        {
            case tools::misc::BLOOM_NONE:      return new BloomNull<T>             ();
            case tools::misc::BLOOM_BASIC:     return new BloomSynchronized<T>     (tai_bloom, nbHash);
            case tools::misc::BLOOM_CACHE:     return new BloomCacheCoherent<T>    (tai_bloom, nbHash);
			case tools::misc::BLOOM_NEIGHBOR:  return new BloomNeighborCoherent<T> (tai_bloom, kmersize, nbHash);
            case tools::misc::BLOOM_DEFAULT:   return new BloomCacheCoherent<T>    (tai_bloom, nbHash);
            default:        throw system::Exception ("bad Bloom kind %d in createBloom", kind);
        }
    }

    /** */
    template<typename T> IBloom<T>* createBloom (
        const std::string& name,
        const std::string& sizeStr,
        const std::string& nbHashStr,
        const std::string& kmerSizeStr
    )
    {
        tools::misc::BloomKind kind;  parse (name, kind);
        return createBloom<T> (kind, (u_int64_t)atol (sizeStr.c_str()), (size_t)atol (nbHashStr.c_str()), atol (kmerSizeStr.c_str()));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_ */
