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
/** \file LargeInt<1>.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */

/** \brief Large integer class
 */





template<>  class LargeInt<1> :  private misc::ArrayData<u_int64_t, 1>
{
public:

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    LargeInt<1> (const u_int64_t& c=0)  {  value[0] = c;  }

     u_int64_t getVal ()  { return *value; }

    static const char* getName ()  { return "LargeInt<1>"; }

    static const size_t getSize ()  { return 8*sizeof(u_int64_t); }

    /** Returns lower 64 bits */
    u_int64_t toInt () const  {  return value[0];  }

    LargeInt<1> operator+  (const LargeInt<1>& other)   const   {  return value[0] + other.value[0];  }
    LargeInt<1> operator-  (const LargeInt<1>& other)   const   {  return value[0] - other.value[0];  }
    LargeInt<1> operator|  (const LargeInt<1>& other)   const   {  return value[0] | other.value[0];  }
    LargeInt<1> operator*  (const int& coeff)           const   {  return value[0] * coeff;        }
    LargeInt<1> operator/  (const u_int32_t& divisor)   const   {  return value[0] / divisor;      }
    u_int32_t   operator%  (const u_int32_t& divisor)   const   {  return value[0] % divisor;      }
    LargeInt<1> operator^  (const LargeInt<1>& other)   const   {  return value[0] ^ other.value[0];  }
    LargeInt<1> operator&  (const LargeInt<1>& other)   const   {  return value[0] & other.value[0];  }
    LargeInt<1> operator&  (const char& other)          const   {  return value[0] & other;        }
    LargeInt<1> operator~  ()                           const   {  return ~value[0];               }
    LargeInt<1> operator<< (const int& coeff)           const   {  return value[0] << coeff;       }
    LargeInt<1> operator>> (const int& coeff)           const   {  return value[0] >> coeff;       }
    bool        operator!= (const LargeInt<1>& c)       const   {  return value[0] != c.value[0];     }
    bool        operator== (const LargeInt<1>& c)       const   {  return value[0] == c.value[0];     }
    bool        operator<  (const LargeInt<1>& c)       const   {  return value[0] < c.value[0];      }
    bool        operator<= (const LargeInt<1>& c)       const   {  return value[0] <= c.value[0];     }

    LargeInt<1>& operator+=  (const LargeInt<1>& other)    {  value[0] += other.value[0]; return *this; }
    LargeInt<1>& operator^=  (const LargeInt<1>& other)    {  value[0] ^= other.value[0]; return *this; }

    u_int8_t  operator[]  (size_t idx) const   {  return (value[0] >> (2*idx)) & 3; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const LargeInt<1> & l)
    {
        s << std::hex << l.value[0] << std::dec;  return s;
    }
    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=32).
     */
    inline void printASCII ( size_t sizeKmer = 32)
    {
        int i;
        u_int64_t temp = value[0];

        
        char seq[33];
        char bin2NT[4] = {'A','C','T','G'};
        
        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }
            seq[sizeKmer]='\0';
        
        std::cout << seq << std::endl;
    }

    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=32).
     */
    std::string toString (size_t sizeKmer) const
    {
        int i;
        u_int64_t temp = value[0];

        char seq[33];
        char bin2NT[4] = {'A','C','T','G'};

        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }
        seq[sizeKmer]='\0';
        return seq;
    }

    
	
	

	
	
    /********************************************************************************/
    inline static u_int64_t revcomp64 (const u_int64_t& x, size_t sizeKmer)
    {
        u_int64_t res = x;

        unsigned char* kmerrev  = (unsigned char *) (&(res));
        unsigned char* kmer     = (unsigned char *) (&(x));

        for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        return (res >> (2*( 32 - sizeKmer))) ;
    }

    /********************************************************************************/
    inline static u_int64_t hash64 (u_int64_t key, u_int64_t seed)
    {
        u_int64_t hash = seed;
        hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
        hash = (~hash) + (hash << 21); // hash = (hash << 21) - hash - 1;
        hash = hash ^ (hash >> 24);
        hash = (hash + (hash << 3)) + (hash << 8); // hash * 265
        hash = hash ^ (hash >> 14);
        hash = (hash + (hash << 2)) + (hash << 4); // hash * 21
        hash = hash ^ (hash >> 28);
        hash = hash + (hash << 31);
		
		return hash;
    }

    /********************************************************************************/
    inline static u_int64_t oahash64 (u_int64_t elem)
    {
        u_int64_t code = elem;
        code = code ^ (code >> 14); //supp
        code = (~code) + (code << 18);
        code = code ^ (code >> 31);
        code = code * 21;
        code = code ^ (code >> 11);
        code = code + (code << 6);
        code = code ^ (code >> 22);
        return code;
    }

	inline static bool comp_less_minimizer (u_int64_t a, u_int64_t b)
	{
		bool res = false;
		
		int mm = system::g_msize;
		u_int64_t * freqm = system::freq_mmer;
		
		//2 revcomp, couteux
		u_int64_t ca = std::min ( NativeInt64::revcomp8(a,mm), a); // revcomp8 version optim revcomp 8nt max
		u_int64_t cb = std::min ( NativeInt64::revcomp8(b,mm), b);
		
		if (freqm[ca] < freqm[cb]   ||
			( (freqm[ca] == freqm[cb]) && (a < b) )
			)
			res = true;
		return res;
		
		//pour minim standard :
		//return (a<b);
	}

	
	/********************************************************************************/
    inline static u_int64_t mhash64 (u_int64_t forw,u_int64_t rev,   u_int64_t prevmin)
    {
		
		u_int64_t res ;
		int ks = system::g_ksize;
		int mm = system::g_msize;
		u_int64_t * freqm = system::freq_mmer;
		u_int64_t maskm = (1<<(2*mm)) - 1;
		
		u_int64_t maskms = (1<<(2*(mm-1))) - 1;
		
		
		int dec = 2*(ks- (mm-1)); // m-1 suffix
		u_int64_t  f_sortant = (forw >> dec ) & maskms; // on a plus le sortant en fait, juste son m-1 suffixe
		u_int64_t  f_new = forw & maskm;
		
		u_int64_t  r_sortant = rev & maskms; //prefix du sortant
		u_int64_t  r_new = (rev >> (2*(ks- (mm))) ) & maskm;
		
		
		// 1000000000 == non init, debut seq
		if(f_sortant== (prevmin & maskms) ||  r_sortant== ((prevmin >>2) & maskms)  || prevmin == 1000000000 ) // peut etre il est sorti; on recalc tout
		{
			u_int64_t minim = 1000000000;
			
			minim = forw  & maskm; // init avec le premier
			for (int i = 0; i <  (1+(ks *2)- 2*mm );i+=2 )
			{
				u_int64_t	 newm = (forw >> i) & maskm;
				if (comp_less_minimizer(newm,minim )  ) minim = newm;
			}
			
			for (int i = 0; i <  (1+(ks *2)- 2*mm );i+=2 )
			{
				u_int64_t	 newm = (rev >> i) & maskm;
				if (comp_less_minimizer(newm,minim )  ) minim = newm;
			}
			
			res = minim;
			
		}
		else //juste besoin comparer aux  2 minim entrant
		{
			res = prevmin;
			
			u_int64_t minim = r_new;
			if (comp_less_minimizer(f_new,r_new )  )
				minim = f_new;
			if (comp_less_minimizer(minim,res )  )
				res = minim;
			
		}
		
		return res;

    }
	
    /********************************************************************************/
    /** computes a simple, naive hash using only 16 bits from input key
     * \param[shift] in : selects which of the input byte will be used for hash computation
     */
    inline static  u_int64_t    simplehash16_64   (u_int64_t key, int  shift)
    {
        u_int64_t input = key >> shift;
        u_int64_t res = random_values[input & 255]   ;
        
        input = input  >> 8;
        res  ^= random_values[input & 255] ;
        
        return res;
        //could be improved by xor'ing result of multiple bytes
    }

    /********************************************************************************/
    static hid_t hdf5 (bool& isCompound)
    {
        return H5Tcopy (H5T_NATIVE_UINT64);
    }

    /********************************************************************************/
    template<typename Map>
    static LargeInt<1> polynom (const char* data, size_t size, Map fct)
    {
        LargeInt<1> res (0);
        for (size_t i=0; i<size; ++i)  {  res.value[0] = 4 * res.value[0] + fct(data[i]);  }
        return res;
    }

private:

    friend LargeInt<1> revcomp (const LargeInt<1>& i,   size_t sizeKmer);
    friend u_int64_t    hash1    (const LargeInt<1>& key, u_int64_t  seed);
    friend u_int64_t    oahash  (const LargeInt<1>& key);
	friend u_int64_t    mhash  (const LargeInt<1>& key,const LargeInt<1>& key2,   u_int64_t prevmin);
    friend u_int64_t    simplehash16    (const LargeInt<1>& key, int  shift);

};

/********************************************************************************/
inline LargeInt<1> revcomp (const LargeInt<1>& x, size_t sizeKmer)
{
    return LargeInt<1>::revcomp64 (x.value[0], sizeKmer);
}

/********************************************************************************/
inline u_int64_t hash1 (const LargeInt<1>& key, u_int64_t seed=0)
{

    return LargeInt<1>::hash64 (key.value[0], seed);
}

/********************************************************************************/
inline u_int64_t oahash (const LargeInt<1>& key)
{
    return LargeInt<1>::oahash64 (key.value[0]);
}

/********************************************************************************/
inline u_int64_t mhash (const LargeInt<1>& key, const LargeInt<1>& key2, u_int64_t prevmin)
{
	//printf("call mash L1  .. ");
    return LargeInt<1>::mhash64 (key.value[0],key2.value[0],prevmin);
}

/********************************************************************************/
inline u_int64_t simplehash16 (const LargeInt<1>& key, int  shift)
{
    return LargeInt<1>::simplehash16_64 (key.value[0], shift);
}
