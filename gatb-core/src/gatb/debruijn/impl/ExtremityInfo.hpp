#ifndef _GATB_CORE_DEBRUIJN_IMPL_EXTREMITY_INFO_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_EXTREMITY_INFO_HPP_
 
/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

    enum Unitig_pos 
    {
        // cannot start at 0 because i'm using (pos & UNITIG_BEGIN) as a test in many places.
        UNITIG_BEGIN = 1,
        UNITIG_END = 2,
        UNITIG_BOTH = 3, /* not encoded in ExtremityInfo, but encoded in NodeGU */
        UNITIG_INSIDE = 4 /* not encoded in ExtremityInfo, but encoded in NodeGU */
    };


   
// core unitigs graph part
    // btw
    // bit packing so dirty i'd need to call O2 to get it cleaned up
    struct ExtremityInfo 
    {
        public:
        uint64_t unitig;
        bool rc; /* when encoding an extremity:
                     whether the kmer in canonical form appears as rc in the unitig
                    when encoding a link:
                     whether you need to reverse the next unitig in the link
                 */
        Unitig_pos pos; // whether the kmer is at left extremity of unitig or right extremity
        ExtremityInfo(uint64_t u, bool r, Unitig_pos p) : unitig(u),rc(r), pos(p) {}
        ExtremityInfo() {} // because i defined another constructor
        ExtremityInfo(const uint64_t val) { unpack(val); } 
        std::string toString() const
        { return " unitig: " + std::to_string(unitig) + " rc:" + std::to_string(rc) + " p:" + ((pos&UNITIG_BEGIN)?"UNITIG_BEGIN":"") + ((pos&UNITIG_END)?"UNITIG_END":"") ; }
        uint64_t pack()
        {
            if ((int)pos > 2) { std::cout << "incorrect encoding for pos in packed ExtremityInfo: " << (int)pos << std::endl; exit(1); }
            return (pos-1) + (rc << 1) + (unitig << 2);
        }
        uint64_t pack_norc() // possibly used in "merci.cpp"
        {
            return (pos-1) + ( 0 << 1) + (unitig << 2);
        }
        void unpack(uint64_t val)
        {
            
            pos = (val&1)?UNITIG_END:UNITIG_BEGIN;
                                   val >>= 1;
            rc = val&1;            val >>= 1;
            unitig = val;
        }
    } ;
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif

