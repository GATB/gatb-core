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

#include <gatb/kmer/impl/DebloomMinimizerAlgorithm.hpp>

#include <gatb/kmer/impl/PartiInfo.hpp>
//
#include <gatb/tools/storage/impl/StorageTools.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage::impl;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
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
DebloomMinimizerAlgorithm<span>::DebloomMinimizerAlgorithm (
    Storage& storage,
    Storage& storageSolids,
    Partition<Count>*    solidIterable,
    size_t              kmerSize,
    size_t              max_memory,
    size_t              nb_cores,
    BloomKind           bloomKind,
    DebloomKind         cascadingKind,
    const std::string&  debloomUri,
    IProperties*        options
)
    :  DebloomAlgorithm<span>(storage, storageSolids, solidIterable, kmerSize, max_memory, nb_cores, bloomKind, cascadingKind, debloomUri, options)
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
DebloomMinimizerAlgorithm<span>::DebloomMinimizerAlgorithm (tools::storage::impl::Storage& storage)
:  DebloomAlgorithm<span>(storage)
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
template<typename Model, typename ModelMini, typename Count, typename Type>
struct FunctorKmersExtension
{
    struct FunctorNeighbors
    {
        IBloom<Type>*         bloom;
        ModelMini&           _modelMini;
        vector<Type>&        _solids;
        PartitionCache<Type> _partition;
        Repartitor&          _repart;

        FunctorNeighbors (
            IBloom<Type>*       bloom,
            ModelMini&          modelMini,
            vector<Type>&       solids,
            Partition<Type>*    extentParts,
            Repartitor&         repart
        )
            : bloom(bloom), _modelMini(modelMini), _solids(solids), _partition(*extentParts,1<<12), _repart(repart)  {}

        void operator() (const Type& neighbor)//  const
        {
            /** We can already get rid of neighbors that are in the current solid kmers partition. */
            if (std::binary_search (_solids.begin(), _solids.end(), neighbor)==false)
            {
                /** We get the partition index of the neighbor from its minimizer value and the
                 * minimizer repartition table.
                 * Note : we have here to compute minimizers from scratch, which may be time expensive;
                 * maybe a better way could be found. */
                u_int64_t mm = _repart (_modelMini.getMinimizerValue(neighbor));
                _partition[mm].insert (neighbor);
            }
        }

    } functorNeighbors;

    Model&        model;
    IBloom<Type>* bloom;

    FunctorKmersExtension (
        Model&              model,
        ModelMini&          modelMini,
        IBloom<Type>*       bloom,
        Partition<Type>*    extentParts,
        vector<Type>&       solids,
        Repartitor&         repart,
        size_t              currentPart
    )
        : model(model),  bloom(bloom),  functorNeighbors(bloom,modelMini, solids, extentParts, repart) {}

    void operator() (const Count& kmer) const
    {
        /** We want to know which neighbors of the current kmer are in the Bloom filter.
         * Note that, according to the Bloom filter implementation, we can have optimized
         * way to get 8 answers in one shot. */
        bitset<8> mask =  bloom->contains8 (kmer.value);

        /** We iterate the neighbors (only those found in the Bloom filter). */
        model.iterateNeighbors (kmer.value, functorNeighbors, mask);
    }
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
struct FinalizeCmd : public ICommand, public SmartPointer
{
    typedef typename kmer::impl::Kmer<span>::Type   Type;
    typedef typename kmer::impl::Kmer<span>::Count  Count;

    size_t currentIdx;
    Iterable<Count>& solids;
    Iterable<Type>& cfp;
    BagCache<Type> result;

    FinalizeCmd (size_t currentIdx, Collection<Count>& solids, Collection<Type>& cfp, Bag<Type>* bag, ISynchronizer* synchro)
        : currentIdx(currentIdx), solids(solids), cfp(cfp), result(bag,8*1024, synchro)
    {}

    void execute ()
    {
        size_t k=0;

        vector<Type> vecCFP (cfp.getNbItems());

        /** We insert all the cfp items into a vector. */
        Iterator<Type>*  itCFP   = cfp.iterator();   LOCAL(itCFP);
        for (itCFP->first(); !itCFP->isDone(); itCFP->next())  { vecCFP[k++] = itCFP->item(); }

        /** We sort this cfp vector. */
        std::sort (vecCFP.begin(), vecCFP.end());

        /** We need two iterators on the two sets (supposed to be ordered). */
        typename vector<Type>::iterator itVecCFP;
        Iterator<Count>* itKmers = solids.iterator();   LOCAL(itKmers);

        /** We initialize the two iterators. */
        itVecCFP = vecCFP.begin();
        itKmers->first();

        /** We compute the difference between the cfp and the solid kmers.
         *  (see http://www.cplusplus.com/reference/algorithm/set_difference)
         */
        while (itVecCFP != vecCFP.end() && !itKmers->isDone())
        {
            if (*itVecCFP < itKmers->item().value)
            {
                result.insert (*itVecCFP);

                Type tmp = *itVecCFP;  while (++itVecCFP != vecCFP.end() && *(itVecCFP)==tmp) { }
            }
            else if (itKmers->item().value < *itVecCFP)
            {
                itKmers->next();
            }
            else
            {
                Type tmp = *itVecCFP;  while (++itVecCFP != vecCFP.end()  && *(itVecCFP)==tmp) { }
                itKmers->next();
            }
        }

        /** We complete the remaining potential cFP items. */
        for ( ;  itVecCFP != vecCFP.end(); ++itVecCFP)  {  result.insert (*itVecCFP);  }
    }
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
void DebloomMinimizerAlgorithm<span>::execute_aux (
    IProperties* bloomProps,
    IProperties* cfpProps,
    u_int64_t&   totalSizeBloom,
    u_int64_t&   totalSizeCFP
)
{
    /** We need two kmer models : one canonical, the other minimizers. */
    typedef typename Kmer<span>::template ModelMinimizer<Model>  ModelMini;
    Model      model     (this->_kmerSize);
    ModelMini  modelMini (this->_kmerSize, this->_miniSize);

    /** We retrieve the minimizers distribution from the solid kmers storage. */
    Repartitor repart;
    repart.load (this->_storageSolids().getGroup("dsk"));

    /** We create the collection that will hold the critical false positive kmers. */
    Collection<Type>* criticalCollection = new CollectionFile<Type> ("cfp");
    LOCAL (criticalCollection);

    /***************************************************/
    /** We create a bloom and insert solid kmers into. */
    /***************************************************/
    IBloom<Type>* bloom = StorageTools::singleton().loadBloom<Type> (this->_groupBloom, "bloom");
    bloom->use ();
    totalSizeBloom = bloom->getBitSize();

    DEBUG (("DebloomMinimizerAlgorithm<span>::execute_aux  totalSizeBloom=%lld \n", totalSizeBloom));

    /** We get the number of partitions in the solid kmers set. */
    size_t nbPartitions = this->_solidIterable->size();

    /** We use a temporary partition that will hold the neighbors extension of the solid kmers. */
    Storage* cfpPartitions = StorageFactory(STORAGE_HDF5).create ("debloom_partitions", true, false);
    LOCAL (cfpPartitions);
    Partition<Type>* debloomParts = & (*cfpPartitions)().getPartition<Type> ("parts", nbPartitions);

    /*************************************************/
    /** We build the solid neighbors extension.      */
    /*************************************************/
    {
        TIME_INFO (this->getTimeInfo(), "fill_debloom_file");

        DEBUG (("DebloomMinimizerAlgorithm<span>::execute_aux   fill_debloom_file BEGIN   nbParts=%ld\n", nbPartitions));

        /** We create an iterator for progress information. */
        Iterator<int>* itParts = this->createIterator (
            new Range<int>::Iterator (0,nbPartitions-1), nbPartitions, DebloomAlgorithm<span>::progressFormat2()
        );
        LOCAL (itParts);

        /** We iterate the [kmer,count] objects partitions*/
        for (itParts->first (); !itParts->isDone(); itParts->next())
        {
            /** Shortcut. */
            size_t p = itParts->item();

            /** We retrieve an iterator on the Count objects of the pth partition. */
            Iterator<Count>* itKmers = (*this->_solidIterable)[p].iterator();
            LOCAL (itKmers);

            /** We fill a vector with kmers only (don't care about counts here).
             * The items in the partition are supposed to be sorted, so will be this vector.
             * THIS IS IMPORTANT BECAUSE we will use a binary search on that vector. */
            vector<Type> solids ( ((*this->_solidIterable)[p]).getNbItems());
            size_t k=0;  for (itKmers->first(); !itKmers->isDone(); itKmers->next()) { solids[k++] = itKmers->item().value; }

            /** We create functor that computes the neighbors extension of the solid kmers. */
            FunctorKmersExtension<Model,ModelMini,Count,Type> functorKmers (model, modelMini, bloom, debloomParts, solids, repart, p);

            /** We iterate the solid kmers. */
            this->getDispatcher()->iterate (itKmers, functorKmers);

        }  /* for (itParts->first (); ...) */

        /** We flush the built partition. */
        debloomParts->flush();

        DEBUG (("DebloomMinimizerAlgorithm<span>::execute_aux   fill_debloom_file END \n"));
    }

    /** We get rid of the bloom. */
    bloom->forget ();

    /*************************************************************/
    /** We extract the solid kmers from the neighbors extension. */
    /*************************************************************/
    {
        TIME_INFO (this->getTimeInfo(), "finalize_debloom_file");

        DEBUG (("DebloomMinimizerAlgorithm<span>::execute_aux   finalize_debloom_file BEGIN \n"));

        ISynchronizer* synchro = System::thread().newSynchronizer();  LOCAL (synchro);

        /** We create an iterator for progress information. */
        Iterator<int>* itParts = this->createIterator (new Range<int>::Iterator (0,nbPartitions-1), nbPartitions, DebloomAlgorithm<span>::progressFormat3());
        LOCAL (itParts);

        size_t nbCoresMax = this->getDispatcher()->getExecutionUnitsNumber();

        for (itParts->first (); !itParts->isDone(); )
        {
            vector<ICommand*> cmd;

            size_t cfpSize = 0;
            size_t nbCores = 0;

            for ( ;  !itParts->isDone() && (nbCores<nbCoresMax && cfpSize < this->_max_memory*MBYTE);  itParts->next())
            {
                /** Shortcut. */
                size_t p = itParts->item();

                nbCores ++;
                cfpSize += (*debloomParts)[p].getNbItems() * sizeof(Type);

                cmd.push_back (new FinalizeCmd<span> (p, (*this->_solidIterable)[p], (*debloomParts)[p], criticalCollection, synchro));
            }

            this->getDispatcher()->dispatchCommands (cmd);
        }

        /** We flush the critical false positive collection. */
        criticalCollection->flush();

        /** We can remove the temporary partition. */
        cfpPartitions->remove();

        DEBUG (("DebloomMinimizerAlgorithm<span>::execute_aux   finalize_debloom_file END \n"));
    }

    /*************************************************************/
    /** We build the final cFP container.                        */
    /*************************************************************/
    {
        TIME_INFO (this->getTimeInfo(), "cascading");

        /** We build the container node. */
        this->createCFP (criticalCollection, cfpProps, totalSizeCFP);
    }

    /** We remove the two temporary files. */
    criticalCollection->remove ();
}

/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomMinimizerAlgorithm <KSIZE_1>;
template class DebloomMinimizerAlgorithm <KSIZE_2>;
template class DebloomMinimizerAlgorithm <KSIZE_3>;
template class DebloomMinimizerAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
