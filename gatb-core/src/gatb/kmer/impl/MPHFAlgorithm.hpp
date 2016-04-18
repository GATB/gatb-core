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

#ifndef _GATB_CORE_KMER_IMPL_MPHF_ALGORITHM_HPP_
#define _GATB_CORE_KMER_IMPL_MPHF_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/api/Enums.hpp> // for MPHFKind
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/MapMPHF.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** \brief Algorithm that builds a hash table whose keys are kmers and values are
 *  kmer abundances.
 *
 *  This class uses a [kmer,abundance] mapping by using a minimal perfect hash function (MPHF).
 *  For N kmers (ie. the keys), the hash function gives a unique integer value between
 *  0 and N-1.
 *
 *  It uses two template parameters:
 *      1) span : gives the max usable size for kmers
 *      2) Abundance_t : type of the abundance values (on 1 byte by default)
 *      2) NodeState_t : type of the node states values (on half a byte by default, grouped by two per byte)
 *
 *  Storing the values (ie. the abundances) is done by creating a vector of size N. Asking
 *  the abundance of a kmer consists in:
 *      * getting the hash code H of the kmer
 *      * getting the object at index H in the vector of values.
 *
 * The MPHF function is built from a list of kmers values of type Kmer<span>::Type. Since
 * the building of the MPHF may take a while, it is saved in a Storage object; more precisely,
 * it is saved in a collection given by a couple [group,name]. Such a couple is likely to be
 * the group of the SortingCount algorithm, with a name being by convention "mphf".
 *
 * Once the MPHF is built, it is populated by the kmers abundance values, which means that we
 * set each value of each key of the hash table. The abundances are clipped to a maximum value
 * in order not to exceed the Abundance_t type capacity (provided as a template of the
 * MPHFAlgorithm class). The maximum value is computed through the std::numeric_limits traits.
 *
 * Once the abundance map is built and populated, it is available through the 'getAbundanceMap' method. It may
 * be used for instance by the Graph class in order to get the abundance of any node (ie. kmer)
 * of the de Bruijn graph.
 *
 * Note: the keys of the hash table are of type Kmer<span>::Type, but we need however to have
 * the abundance information through the Kmer<span>::Count type. That's why we need to use
 * 2 Iterable instances, one of type Kmer<span>::Count and one of type Kmer<span>::Type.
 *
 * Some statistics about the MPHF building are gathered and put into the Properties 'info'.
 */
template<size_t span=KMER_DEFAULT_SPAN, typename Abundance_t=u_int8_t, typename NodeState_t=u_int8_t>
class MPHFAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::Type  Type;
    typedef typename kmer::impl::Kmer<span>::Count Count;

    /** We define the maximum abundance according to the provided type (value set in the cpp file). */
    static const Abundance_t MAX_ABUNDANCE;

    /** We define the type of the hash table of couples [kmer/abundance]. */
    typedef tools::collections::impl::MapMPHF<Type,Abundance_t>  AbundanceMap;
    
    /** We define the type of the hash table of couples [kmer/node state]. */
    typedef tools::collections::impl::MapMPHF<Type,NodeState_t>  NodeStateMap;

    /** We define the type of the hash table of couples [kmer/graph adjacency information]. */
    typedef u_int8_t Adjacency_t;
    typedef tools::collections::impl::MapMPHF<Type,Adjacency_t>  AdjacencyMap;


    /** Constructor.
     * \param[in] group : storage group where to save the MPHF once built
     * \param[in] name : name of the collection in the group where the MPHF will be saved
     * \param[in] solidCounts : iterable on couples [kmers/abundance]
     * \param[in] solidKmers  : iterable on kmers
     * \param[in] buildOrLoad : true for build/save the MPHF, false for load only
     * \param[in] options : extra options for configuration (may be empty) */
    MPHFAlgorithm (
        tools::misc::MPHFKind                 mphfKind,
        tools::storage::impl::Group&          group,
        const std::string&                    name,
        tools::collections::Iterable<Count>*  solidCounts,
        tools::collections::Iterable<Type>*   solidKmers,
        unsigned int                          nbCores,
        bool                                  buildOrLoad,
        tools::misc::IProperties*   options    = 0
    );

    /** Destructor. */
    ~MPHFAlgorithm ();

    /** Implementation of the Algorithm::execute method. */
    void execute ();

    /** Get the number of bits of a value.
     * \return the number of bits per kmer. */
    float getNbBitsPerKmer () const;

    /** Accessor to the map. Note : if clients get this map and use it (as a SmartPointer),
     * the map instance will be still alive (ie. not deleted) even if the MPHFAlgorithm
     * instance that built it is deleted first.
     * \return the map instance. */
    AbundanceMap* getAbundanceMap () const  { return _abundanceMap; }
    NodeStateMap* getNodeStateMap () const  { return _nodeStateMap; }
    NodeStateMap* getAdjacencyMap () const  { return _adjacencyMap; }

private:

    tools::storage::impl::Group& _group;
    std::string                  _name;
    bool                         _buildOrLoad;
    size_t                       _dataSize;
    size_t                       _nb_abundances_above_precision;

    /** Iterable on the couples [kmer,abundance] */
    tools::collections::Iterable<Count>* _solidCounts;
    void setSolidCounts (tools::collections::Iterable<Count>* solidCounts)  {  SP_SETATTR(solidCounts); }

    /** Iterable on the kmers */
    tools::collections::Iterable<Type>* _solidKmers;
    void setSolidKmers (tools::collections::Iterable<Type>* solidKmers)  {  SP_SETATTR(solidKmers); }

    /** Hash table instance. */
    AbundanceMap* _abundanceMap;
    NodeStateMap* _nodeStateMap;
    AdjacencyMap* _adjacencyMap;
    void setAbundanceMap (AbundanceMap* abundanceMap)  { SP_SETATTR(abundanceMap); }
    void setNodeStateMap (NodeStateMap* nodeStateMap)  { SP_SETATTR(nodeStateMap); }
    void setAdjacencyMap (AdjacencyMap* adjacencyMap)  { SP_SETATTR(adjacencyMap); }

    /** Set the abundance for each entry in the hash table. */
    void populate ();
    
    /** Initialize the node state for each entry in the hash table. */
    void initNodeStates ();

    /** Check the content of the map once built. */
    void check ();

    /** We define a specific Progress class for progress feedback during hash function building.
     * We need a special implementation here because of emphf (we can't modify too much code in
     * emphf, so we have to hack it some stuff here). */
    class ProgressCustom : public tools::misc::impl::ProgressProxy
    {
    public:
        ProgressCustom (tools::dp::IteratorListener* ref);  void reset (u_int64_t ntasks);
    private:
        size_t nbReset;
    };

    /** Progress instance. */
    tools::dp::IteratorListener* _progress;
    void setProgress (tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }
	
	/**  rememer mphf kind here also*/
	tools::misc::MPHFKind _mphfKind;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MPHF_ALGORITHM_HPP_ */

