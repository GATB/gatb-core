#ifdef WITH_MPHF
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

#ifndef _MPHF_HPP_
#define _MPHF_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <gatb/debruijn/api/IContainerNode.hpp>

#include <gatb/tools/collections/impl/MPHF.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class MPHFAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:
    typedef unsigned char mphf_abundance_t;

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::Type  Type;
    typedef typename kmer::impl::Kmer<span>::Count Count;
    typedef typename gatb::core::tools::collections::impl::MPHF<span, mphf_abundance_t> MPHF;

    /** */
    MPHFAlgorithm (
        tools::storage::impl::Storage& storage,
        tools::collections::Iterable<Count>* solidIterable,
        size_t                      kmerSize,
        const std::string&          mphfUri = "mphf",
        tools::misc::IProperties*   options    = 0,
        bool construct_emphf = false 
    );

    /** */
    ~MPHFAlgorithm ();

    /** */
    void execute ();

    /** */
    float getNbBitsPerKmer () const;


    MPHF *getMPHF();


private:

    /** */
    tools::storage::impl::Storage& _storage;

    tools::storage::impl::Group& _group;

    size_t       _kmerSize;

    std::string  _mphfUri;

    tools::collections::Iterable<Count>* _solidIterable;
    void setSolidIterable (tools::collections::Iterable<Count>* solidIterable)  {  SP_SETATTR(solidIterable); }

    MPHF *mphf_class;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
#endif
#endif 

