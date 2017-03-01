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

#include <gatb/debruijn/impl/UnitigsConstructionAlgorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/bcalm2/bcalm_algo.hpp>
#include <gatb/bcalm2/bglue_algo.hpp>

#include <queue>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace debruijn  {   namespace impl {
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
UnitigsConstructionAlgorithm<span>::UnitigsConstructionAlgorithm (
    tools::storage::impl::Storage& storage,
    std::string                 unitigs_filename,
    size_t                      nb_cores,
    tools::misc::IProperties*   options
)
    : Algorithm("bcalm2-wrapper", nb_cores, options), _storage(storage), unitigs_filename(unitigs_filename)
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
UnitigsConstructionAlgorithm<span>::~UnitigsConstructionAlgorithm ()
{
}

/*********************************************************************/
/*********************************************************************/

template <size_t span>
void UnitigsConstructionAlgorithm<span>::execute ()
{
    int kmerSize =
            getInput()->getInt(STR_KMER_SIZE);
    int abundance = 
            getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
    int minimizerSize =
        getInput()->getInt(STR_MINIMIZER_SIZE);
    int nb_threads =
        getInput()->getInt(STR_NB_CORES);
    int minimizer_type =
        getInput()->getInt(STR_MINIMIZER_TYPE);
    bool verbose = getInput()->getInt(STR_VERBOSE);
        
    unsigned int nbThreads = this->getDispatcher()->getExecutionUnitsNumber();
    if ((unsigned int)nb_threads > nbThreads)
    {
        std::cout << "Uh. Unitigs graph construction called with nb_threads " << nb_threads << " but dispatcher has nbThreads " << nbThreads << std::endl;

    }

    bcalm2<span>(&_storage, unitigs_filename, kmerSize, abundance, minimizerSize, nbThreads, minimizer_type, verbose); 
    bglue<span> (&_storage, unitigs_filename, kmerSize,            minimizerSize, nbThreads, minimizer_type, verbose);
    link_unitigs(unitigs_filename, kmerSize, verbose);

    /** We gather some statistics. */
    //getInfo()->add (1, "stats");
    //getInfo()->add (2, "nb_unitigs", "%ld", /* */);
    
    //getInfo()->add (1, "time");
    //getInfo()->add (2, "build", "%.3f", /* */);
}

/* this procedure finds the overlaps between unitigs, using a hash table of all extremity (k-1)-mers
 *
 * I guess it's like AdjList in ABySS. It's also like contigs_to_fastg in MEGAHIT.
 * 
 * could be replaced by keeping edges during BCALM2, but it's not the case for now */
template<size_t span>
void UnitigsConstructionAlgorithm<span>::
link_unitigs(string unitigs_filename, int kmerSize, bool verbose)
{
    bool debug = false;

    BankFasta inputBank (unitigs_filename);
    BankFasta::Iterator itSeq (inputBank);
    uint64_t utig_counter = 0;
    
    Model modelKminusOne(kmerSize - 1); // it's canonical (defined in the .hpp file)

    if (verbose)
        std::cout << "Finding links between unitigs, pass 1,  mem current/maxRSS: " << system::impl::System::info().getMemorySelfUsed() / 1024 << "/" << system::impl::System::info().getMemorySelfMaxUsed() / 1024  << std::endl;

    for (itSeq.first(); !itSeq.isDone(); itSeq.next()) 
    {
        const string& seq = itSeq->toString();
 
        typename Model::Kmer kmerBegin = modelKminusOne.codeSeed(seq.substr(0, kmerSize-1).c_str(), Data::ASCII);
        typename Model::Kmer kmerEnd = modelKminusOne.codeSeed(seq.substr(seq.size() - kmerSize+1).c_str(), Data::ASCII);

        bool beginInSameOrientation = modelKminusOne.toString(kmerBegin.value()) == seq.substr(0,kmerSize-1);
        bool endInSameOrientation = modelKminusOne.toString(kmerEnd.value()) == seq.substr(seq.size() - kmerSize+1);

        ExtremityInfo eBegin(utig_counter, !beginInSameOrientation /* because we record rc*/, UNITIG_BEGIN);
        ExtremityInfo eEnd(  utig_counter, !endInSameOrientation,                             UNITIG_END);
        // there is no UNITIG_BOTH here because we're taking (k-1)-mers.

        utigs_links_map[kmerBegin.value()].push_back(eBegin.pack());
        utigs_links_map[kmerEnd.value()].push_back(eEnd.pack());
        utig_counter++;
    }

    BankFasta* out = new BankFasta(unitigs_filename+".indexed");
    
    if (verbose)
        std::cout << "Finding links between unitigs, pass 2,  mem current/maxRSS: " << system::impl::System::info().getMemorySelfUsed() / 1024 << "/" << system::impl::System::info().getMemorySelfMaxUsed() / 1024  << std::endl;

    uint64_t utigs_number = 0;
    for (itSeq.first(); !itSeq.isDone(); itSeq.next()) 
    {
        const string& seq = itSeq->toString();
        const string& comment = itSeq->getComment();
 
        typename Model::Kmer kmerBegin = modelKminusOne.codeSeed(seq.substr(0, kmerSize-1).c_str(), Data::ASCII);
        typename Model::Kmer kmerEnd = modelKminusOne.codeSeed(seq.substr(seq.size() - kmerSize+1).c_str(), Data::ASCII);
        bool beginInSameOrientation =  modelKminusOne.toString(kmerBegin.value()) == seq.substr(0,kmerSize-1); // that could be optimized, revcomp was already computed during codeSeed
        bool endInSameOrientation =  modelKminusOne.toString(kmerEnd.value()) == seq.substr(seq.size() - kmerSize+1);

        // treat special palindromic kmer cases
        bool nevermindInOrientation = false;
        bool nevermindOutOrientation = false;
        if (((kmerSize - 1) % 2) == 0) 
        {
            if (kmerBegin.isPalindrome()) nevermindInOrientation = true;
            if (kmerEnd.isPalindrome())   nevermindOutOrientation = true;

        }
        
        if (debug) std::cout << "unitig: " << seq << std::endl;

        string links;

        // in-neighbors
        for (auto in_packed : utigs_links_map[kmerBegin.value()])
        {
            ExtremityInfo e_in(in_packed);


            if (debug) std::cout << "extremity " << modelKminusOne.toString(kmerBegin.value()) << " ";
            if (debug) std::cout << "potential in-neighbor: " << e_in.toString() << " beginSameOrientation " << beginInSameOrientation;

            // what we want are these four cases:
            //  ------[end same orientation] -> [begin same orientation]----
            //  [begin diff orientation]---- -> [begin same orientation]----
            //  ------[end diff orientation] -> [begin diff orientation]----
            //  [begin same orientation]---- -> [begin diff orientation]----
            if ((((beginInSameOrientation)  &&  (e_in.pos == UNITIG_END  ) && (e_in.rc == false)) ||
                ((beginInSameOrientation) &&  (e_in.pos == UNITIG_BEGIN) && (e_in.rc == true)) ||
              (((!beginInSameOrientation)) && (e_in.pos == UNITIG_END  ) && (e_in.rc == true)) ||
              (((!beginInSameOrientation)) && (e_in.pos == UNITIG_BEGIN) && (e_in.rc == false)))
                    || nevermindInOrientation)
            {
                if (nevermindInOrientation && (e_in.unitig == utigs_number)) continue; // don't consider the same extremity

                //LinkInfo li(e_in.unitig, e_in.rc ^ beginInSameOrientation);
                //incoming[utig_number].push_back(li.pack());
                bool rc = e_in.rc ^ (!beginInSameOrientation);
                links += "L:-:" + to_string(e_in.unitig) + ":" + (rc?"+":"-") + " "; /* invert-reverse because of incoming orientation. it's very subtle and i'm still not sure i got it right */

                if (nevermindInOrientation)
                    links += "L:-:" + to_string(e_in.unitig) + ":" + ((!rc)?"+":"-") + " "; /* in that case, there is also another link with the reverse direction*/
                    
                if (debug) std::cout << " [valid] ";
            }
                
            if (debug) std::cout << std::endl;
        }

        // out-neighbors
        for (auto out_packed : utigs_links_map[kmerEnd.value()])
        {
            ExtremityInfo e_out(out_packed);

            if (debug) std::cout << "extremity " << modelKminusOne.toString(kmerEnd.value()) << " ";
            if (debug) std::cout << "potential out-neighbor: " << e_out.toString();

            // what we want are these four cases:
            //  ------[end same orientation] -> [begin same orientation]----
            //  ------[end same orientation] -> ------[end diff orientation]
            //  ------[end diff orientation] -> [begin diff orientation]----
            //  ------[end diff orientation] -> ------[end same orientation]
            if ((((endInSameOrientation) && (e_out.pos == UNITIG_BEGIN) && (e_out.rc == false)) ||
                ((endInSameOrientation) && (e_out.pos == UNITIG_END  ) && (e_out.rc == true)) ||
             (((!endInSameOrientation)) && (e_out.pos == UNITIG_BEGIN) && (e_out.rc == true)) ||
             (((!endInSameOrientation)) && (e_out.pos == UNITIG_END  ) && (e_out.rc == false)))
                ||nevermindOutOrientation)
            {
                if (nevermindOutOrientation && (e_out.unitig == utigs_number)) continue; // don't consider the same extremity

                //LinkInfo li(e_out.unitig, e_out.rc ^ endInSameOrientation);
                //outcoming[utig_number].push_back(li.pack());
                bool rc = e_out.rc ^ (!endInSameOrientation);
                links += "L:+:" + to_string(e_out.unitig) + ":" + (rc?"-":"+") + " "; /* logically this is going to be opposite of the line above */

                if (nevermindOutOrientation)
                    links += "L:+:" + to_string(e_out.unitig) + ":" + ((!rc)?"-":"+") + " "; /* in that case, there is also another link with the reverse direction*/

                if (debug) std::cout << " [valid] ";
            }
            if (debug) std::cout << std::endl;
        }

        Sequence s (Data::ASCII);
        s.getData().setRef ((char*)seq.c_str(), seq.size());
        s._comment = comment + " " + links;
        out->insert(s);
        utigs_number++;
    }
    nb_unitigs = utigs_number; 

    delete out;
    system::impl::System::file().remove (unitigs_filename);
    system::impl::System::file().rename (unitigs_filename+".indexed", unitigs_filename);
}




/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
