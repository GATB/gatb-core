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

/** \file LargeInt<2>.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */

/********************************************************************************/
#if  INT128_FOUND == 1
/********************************************************************************/

template<>  class LargeInt<2> :  private misc::ArrayData<__uint128_t, 1>
{
public:

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    LargeInt<2> (const __uint128_t& c=0)  {  value[0] = c;  }

     u_int64_t getVal () const  { return *value; }

    static const char* getName ()  { return "LargeInt<2>"; }

    static const size_t getSize ()  { return 8*sizeof(__uint128_t); }

    LargeInt<2> operator+  (const LargeInt<2>& other)     const   {  return value[0] + other.value[0];  }
    LargeInt<2> operator-  (const LargeInt<2>& other)     const   {  return value[0] - other.value[0];  }
    LargeInt<2> operator|  (const LargeInt<2>& other)     const   {  return value[0] | other.value[0];  }
    LargeInt<2> operator*  (const int& coeff)              const   {  return value[0] * coeff;        }
    LargeInt<2> operator/  (const u_int32_t& divisor)      const   {  return value[0] / divisor;      }
    u_int32_t    operator%  (const u_int32_t& divisor)      const   {  return value[0] % divisor;      }
    LargeInt<2> operator^  (const LargeInt<2>& other)     const   {  return value[0] ^ other.value[0];  }
    LargeInt<2> operator&  (const LargeInt<2>& other)     const   {  return value[0] & other.value[0];  }
    LargeInt<2> operator&  (const char& other)             const   {  return value[0] & other;        }
    LargeInt<2> operator~  ()                              const   {  return ~value[0];               }
    LargeInt<2> operator<< (const int& coeff)              const   {  return value[0] << coeff;       }
    LargeInt<2> operator>> (const int& coeff)              const   {  return value[0] >> coeff;       }
    bool         operator!= (const LargeInt<2>& c)         const   {  return value[0] != c.value[0];     }
    bool         operator== (const LargeInt<2>& c)         const   {  return value[0] == c.value[0];     }
    bool         operator<  (const LargeInt<2>& c)         const   {  return value[0] < c.value[0];      }
    bool         operator<= (const LargeInt<2>& c)         const   {  return value[0] <= c.value[0];     }

    LargeInt<2>& operator+=  (const LargeInt<2>& other)    {  value[0] += other.value[0]; return *this; }
    LargeInt<2>& operator^=  (const LargeInt<2>& other)    {  value[0] ^= other.value[0]; return *this; }

    LargeInt<2>& operator<<=  (const int& coeff)  { value[0] <<= coeff; return *this; } 
    LargeInt<2>& operator>>=  (const int& coeff)  { value[0] >>= coeff; return *this; }

    u_int8_t  operator[]  (size_t idx) const   {  return (value[0] >> (2*idx)) & 3; }

    /** Output stream overload. NOTE: for easier process, dump the value in hexadecimal.
     * \param[in] os : the output stream
     * \param[in] in : the integer value to be output.
     * \return the output stream.
     */
    friend std::ostream & operator<<(std::ostream & os, const LargeInt<2> & in)
    {
        __uint128_t x = in.value[0];

        u_int64_t high_nucl = (u_int64_t) (x>>64);
        u_int64_t low_nucl  = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));

        if (high_nucl == 0) {   os << std::hex <<                     low_nucl << std::dec;  }
        else                {   os << std::hex << high_nucl << "." << low_nucl << std::dec;  }
        return os;
    }

    /********************************************************************************/
    
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=64).
     */
    inline void printASCII ( size_t sizeKmer = 64)
    {
        int i;
        u_int64_t temp = value[0];
        
        
        char seq[65];
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
        char seq[65];
        char bin2NT[4] = {'A','C','T','G'};

        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }
        seq[sizeKmer]='\0';
        return seq;
    }

    /********************************************************************************/
    static hid_t hdf5 (bool& isCompound)
    {
        hid_t result = H5Tcopy (H5T_NATIVE_INT);
        H5Tset_precision (result, 128);
        return result;
    }
    
    /********************************************************************************/
    template<typename Map>
    static LargeInt<2> polynom (const char* data, size_t size, Map fct)
    {
        LargeInt<2> res (0);
        for (size_t i=0; i<size; ++i)  {  res.value[0] = res.value[0] * 4 + fct(data[i]);  }
        return res;
    }

private:
    friend LargeInt<2> revcomp (const LargeInt<2>& i,   size_t sizeKmer);
    friend u_int64_t    hash1    (const LargeInt<2>& key, u_int64_t  seed);
    friend u_int64_t    oahash  (const LargeInt<2>& key);
    friend u_int64_t    simplehash16    (const LargeInt<2>& key, int  shift);

};

/********************************************************************************/
inline LargeInt<2> revcomp (const LargeInt<2>& in, size_t sizeKmer)
{
    //                  ---64bits--   ---64bits--
    // original kmer: [__high_nucl__|__low_nucl___]
    //
    // ex:            [         AC  | .......TG   ]
    //
    //revcomp:        [         CA  | .......GT   ]
    //                 \_low_nucl__/\high_nucl/

    const __uint128_t& x = in.value[0];

    u_int64_t high_nucl = (u_int64_t)(x>>64);
    int nb_high_nucl = sizeKmer>32?sizeKmer - 32:0;

    __uint128_t revcomp_high_nucl = NativeInt64::revcomp64 (high_nucl, nb_high_nucl);

    if (sizeKmer<=32) revcomp_high_nucl = 0; // srsly dunno why this is needed. gcc bug? u_int64_t x ---> (x>>64) != 0

    u_int64_t low_nucl = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));
    int nb_low_nucl = sizeKmer>32?32:sizeKmer;

    __uint128_t revcomp_low_nucl = NativeInt64::revcomp64 (low_nucl, nb_low_nucl);

    return (revcomp_low_nucl<<(2*nb_high_nucl)) + revcomp_high_nucl;
}

/********************************************************************************/
inline u_int64_t hash1 (const LargeInt<2>& item, u_int64_t seed=0)
{
    const __uint128_t& elem = item.value[0];

    return NativeInt64::hash64 ((u_int64_t)(elem>>64),seed) ^
           NativeInt64::hash64 ((u_int64_t)(elem&((((__uint128_t)1)<<64)-1)),seed);
}

/********************************************************************************/
inline u_int64_t oahash (const LargeInt<2>& item)
{
    const __uint128_t& elem = item.value[0];

    return NativeInt64::oahash64 ((u_int64_t)(elem>>64)) ^
           NativeInt64::oahash64 ((u_int64_t)(elem&((((__uint128_t)1)<<64)-1)));
}



/********************************************************************************/
inline u_int64_t simplehash16 (const LargeInt<2>& key, int  shift)
{
    return NativeInt64::simplehash16_64 ((u_int64_t)key.value[0], shift);
}

/********************************************************************************/
#endif //INT128_FOUND
/********************************************************************************/
