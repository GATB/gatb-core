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

#ifndef _GATB_CORE_KMER_ICOUNT_PROCESSOR_HPP_
#define _GATB_CORE_KMER_ICOUNT_PROCESSOR_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/Configuration.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
/********************************************************************************/

/** \brief Interface that manages information sent by SortingCountAlgorithm.
 */
template<size_t span>
class ICountProcessor : public system::SmartPointer
{
public:

    /** Shortcuts. */
    typedef typename impl::Kmer<span>::Type Type;

    /** Notification that a [kmer,counts] is available and can be handled by the count processor.
     * \param[in] partId : index of the current partition
     * \param[in] kmer : kmer for which we are receiving counts
     * \param[in] count : vector of counts of the kmer, one count per bank
     * \param[in] sum : sum of the occurrences for all bank.
     */
    virtual bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0) = 0;

    /********************************************************************/
    /*   NOTIFS AT BEGIN/END OF THE MAINLOOP OF SortingCountAlgorithm.  */
    /********************************************************************/

    /** Called just before the mainloop of SortingCountAlgorithm.
     * \param[in] config : configuration of the SortingCountAlgorithm. */
    virtual void begin (const impl::Configuration& config) = 0;

    /** Called just after the mainloop of SortingCountAlgorithm. */
    virtual void end   () = 0;

    /*******************************************************************************/
    /*   NOTIFS AT BEGIN/END OF ONE KMERS PARTITION DURING SortingCountAlgorithm.  */
    /*******************************************************************************/

    /** Clone the instance.
     * An instance can be cloned N times in order to use the cloned instance in one thread.
     * \return the cloned instance. */
    virtual ICountProcessor* clone () = 0;

    /** Called at the beginning of a new kmers partition processing.
     * \param[in] passId : index of the current pass in the SortingCountAlgorithm.
     * \param[in] passId : index of the current kmers partition in the SortingCountAlgorithm.
     * \param[in] cacheSize : memory size used for the current kmers partition
     * \param[in] name : class name of the child PartitionsCommand class.
     */
    virtual void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name) = 0;

    /** Called at the end of a new kmers partition processing.
     * \param[in] passId : index of the current pass in the SortingCountAlgorithm.
     * \param[in] passId : index of the current kmers partition in the SortingCountAlgorithm.
     */
    virtual void endPart   (size_t passId, size_t partId) = 0;

    /*****************************************************************/
    /*                          MISCELLANOUS.                        */
    /*****************************************************************/

    /** Get some properties about the count processor.
     * \return properties. */
    virtual tools::misc::impl::Properties getProperties() const = 0;

    /** Get a vector of instances in case of the current object is a composite.
     * \return a vector of ICountProcessor instance. */
    virtual std::vector<ICountProcessor*> getInstances () const = 0;

    /** Try to get an instance of a specific type within the current object.
     * \return a T pointer to the instance if found, 0 otherwise. */
    template<typename T> T* get () const
    {
        std::vector<ICountProcessor*> v = this->getInstances();
        for (size_t i=0; i<v.size(); i++)  {  if (T* object = dynamic_cast<T*> (v[i]))  { return object; }  }
        return (T*)0;
    }

protected:

    /** Compute the sum of the counts in the provided CountVector. This may be useful to compute
     * this sum only once and provide it several times.
     * \return the sum of occurrences
     */
    virtual CountNumber computeSum (const CountVector& count) const = 0;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_ICOUNT_PROCESSOR_HPP_ */
