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
        UNITIG_BEGIN = 0,
        UNITIG_END = 1,
        UNITIG_BOTH = 2, /* not encoded in ExtremityInfo, but encoded in NodeGU */
        UNITIG_INSIDE = 3 /* not encoded in ExtremityInfo, but encoded in NodeGU */
    };


   
// core unitigs graph part
    // btw
    // bit packing so dirty i'd need to call O2 to get it cleaned up
    struct ExtremityInfo 
    {
        public:
        uint64_t unitig;
        bool deleted;
        bool rc; /* when encoding an extremity:
                     whether the kmer in canonical form appears as rc in the unitig
                    when encoding a link:
                     whether you need to reverse the next unitig in the link
                 */
        Unitig_pos pos; // whether the kmer is at left extremity of unitig or right extremity
        ExtremityInfo(uint64_t u, bool d, bool r, Unitig_pos p) : unitig(u),deleted(d),rc(r), pos(p) {}
        ExtremityInfo() {} // because i defined another constructor
        ExtremityInfo(const uint64_t val) { unpack(val); } 
        std::string toString() const
        { return " rc:" + std::to_string(rc) + " p:" + ((pos&UNITIG_BEGIN)?"left":"") + ((pos&UNITIG_END)?"right":"") + " " + " d:" + std::to_string(deleted); }
        uint64_t pack()
        {
            return pos + (rc << 2) + (deleted << 3) + (unitig << 4);
        }
        void unpack(uint64_t val)
        {
            
            pos = val?UNITIG_END:UNITIG_BEGIN;
                                   val >>= 1;
            rc = val&1;            val >>= 1;
            deleted = val&1;       val >>= 1;
            unitig = val;
        }
    } ;
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif

