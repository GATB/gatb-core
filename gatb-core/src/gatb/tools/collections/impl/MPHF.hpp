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

/** \file MPHF.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Minimal Perfect Hash Function
 */

#ifndef _GATB_CORE_TOOLS_MISC_MPHF_HPP_
#define _GATB_CORE_TOOLS_MISC_MPHF_HPP_

/********************************************************************************/

#include <gatb/system/api/Exception.hpp>
#include <gatb/system/api/types.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

/********************************************************************************/
namespace gatb        {
namespace core        {
namespace tools       {
namespace collections {
namespace impl        {
/********************************************************************************/

typedef std::pair<u_int8_t const*, u_int8_t const*> byte_range_t;

/** For some specialization (see below), we need to adapt the key type to some
 * range of raw data in memory. We provide here a default adaptor that can
 * be used as default template type for the MPHF class.
 */
template<typename T>
struct AdaptatorDefault
{
    byte_range_t operator() (const T& t) const
    {
        const u_int8_t* buf = reinterpret_cast <u_int8_t const*> (&t);
        const u_int8_t* end = buf + sizeof(T);
        return byte_range_t(buf, end);
    }
};

/********************************************************************************/

/** \brief Perfect minimal hash function for a given kind of key
 *
 * This class provides an interface for getting hash codes for some key type T, which
 * can be done through the operator() method
 *
 * This class is not a classic hash feature because it hashes only a given set of T items
 * (provided as a T iterator) through its 'build' method. Once building is done, hash code
 * can be accessed through the operator()
 *
 * We propose here a default implementation that doesn't do much. The idea behind is that
 * we can specialize the class for the 'exist' template argument in order to provide a true
 * implementation (through EMPHF library for instance). If such an implementation exists,
 * the constant 'enabled' will be true, which allows to test it in the code (it is a little
 * bit better than using compilation flag).
 */
template<typename Key, typename Adaptator=AdaptatorDefault<Key>, class Progress=tools::misc::impl::ProgressNone, bool exist=true>
class MPHF : public system::SmartPointer
{
public:
    /** Constant telling whether the feature is enabled or not.
     *
     *   - if not enabled, calls to methods may return an exception.
     *
     *   - if enabled, the implementation should be done through a specialization
     *     of the 'exist' template parameter; such an implementation can be conditionally
     *     compiled through a compilation flag, so if the implementation is not available
     *     on some os/architecture, we switch back to the 'not enabled case' (see EMPHF
     *     for instance)
     */
    static const bool enabled = false;

    /** Definition of a hash value. */
    typedef u_int64_t Code;

    tools::misc::MPHFKind mphfKind;
    
    /** Constructor. */
    MPHF (tools::misc::MPHFKind mphfKind) : mphfKind(mphfKind)   {}
    
    /** Constructor. */
    MPHF () : mphfKind(tools::misc::MPHF_BOOPHF /* boophf is best mphf; so it deserves default*/)   {}

    /** Constructor. */
    MPHF (tools::dp::Iterator<Key>* iterator, size_t nbElts)  { }

    /** Constructor. */
    MPHF (tools::collections::Iterable<Key>* iterable)  { }

    /** Build the hash function from a set of items.
     * \param[in] iterable : keys iterator
     * \param[in] nb threads
     * \param[in] progress : object that listens to the event of the algorithm */
    void build (tools::collections::Iterable<Key>* iterable, int nbThreads = 1, tools::dp::IteratorListener* progress=0)   { error(); }
    
    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Code operator () (const Key& key) { error(); return Code(); }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { error(); return 0; }

    /** Load hash function from a collection
     * \return the number of keys. */
    size_t load (tools::storage::impl::Group& group, const std::string& name) { error(); return 0; }

    /** Save hash function to a collection
     * \return the number of bytes of the saved data. */
    size_t save (tools::storage::impl::Group& group, const std::string& name) { error(); return 0; }

private:

    /** Default error management. */
    void error () const { printf("MPHF error\n"); throw gatb::core::system::ExceptionNotImplemented(); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

/** NOW HERE THE TRICK... We include the wrapper (for EMPHF and BooPHF implementations) if allowed by compilation flag. */
#ifdef WITH_MPHF
    #include <gatb/tools/collections/impl/MPHFWrapper.hpp>
#endif
/* sooo many indirections for the MPHF code.. MPHFAlgorithm > MapMPHF > MPHF > MPHFWrapper > EMPHF/BooMHF. could perhaps simplify someday? */

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_MPHF_HPP_ */
