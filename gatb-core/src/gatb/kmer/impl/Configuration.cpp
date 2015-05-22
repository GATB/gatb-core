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

#include <gatb/kmer/impl/Configuration.hpp>
#include <gatb/system/api/IMemory.hpp>

/********************************************************************************/

using namespace std;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace kmer          {
namespace impl          {
/********************************************************************************/

#define DEBUG(a)  //printf a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Properties Configuration::getProperties() const
{
    Properties result ("config");

    result.add (1, "kmer_size",         "%ld", _kmerSize);
    result.add (1, "mini_size",         "%ld", _minim_size);
    result.add (1, "solidity_kind",      "%s", toString(_solidityKind).c_str());

    std::stringstream ss;
    for (size_t i=0; i<_abundanceUserNb; i++)  {  ss << (i==0 ? "" : " ") << _abundance[i].getBegin(); }

    result.add (1, "abundance_min",     ss.str());
    result.add (1, "abundance_max",     "%ld", _abundance[0].getEnd());

    result.add (1, "available_space",   "%ld", _available_space);
    result.add (1, "sequence_number",   "%ld", _estimateSeqNb);
    result.add (1, "sequence_volume",   "%ld", _estimateSeqTotalSize / system::MBYTE);
    result.add (1, "kmers_number",      "%ld", _kmersNb);
    result.add (1, "kmers_volume",      "%ld", _volume);
    result.add (1, "max_disk_space",    "%ld", _max_disk_space);
    result.add (1, "max_memory",        "%ld", _max_memory);
    result.add (1, "nb_passes",         "%d",  _nb_passes);
    result.add (1, "nb_partitions",     "%d",  _nb_partitions);
    result.add (1, "nb_bits_per_kmer",  "%d",  _nb_bits_per_kmer);
    result.add (1, "nb_cores",          "%d",  _nbCores);
    result.add (1, "partition_type",    "%d",  _partitionType);
    result.add (1, "minimizer_type",    "%s",  (_minimizerType == 0) ? "lexicographic (kmc2 heuristic)" : "frequency");
    result.add (1, "repartition_type",  "%s",  (_repartitionType == 0) ? "unordered" : "ordered");

    result.add (1, "nb_cores_per_partition",     "%d",  _nbCores_per_partition);
    result.add (1, "nb_partitions_in_parallel",  "%d",  _nb_partitions_in_parallel);

    result.add (1, "nb_banks",  "%d",  _nb_banks);

    return result;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

