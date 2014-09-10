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

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

typedef std::pair<uint8_t const*, uint8_t const*> byte_range_t;

/** For some specialization (see below), we need to adapt the key type to some
 * range of raw data in memory. We provide here a default adaptator that can
 * be used as default template type for the MPHF class.
 */
template<typename T>
struct AdaptatorDefault
{
    byte_range_t operator() (const T& t) const
    {
        const uint8_t* buf = reinterpret_cast <uint8_t const*> (&t);
        const uint8_t* end = buf + sizeof(T);
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
template<typename Key, typename Adaptator=AdaptatorDefault<Key>, bool exist=true>
class MPHF
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
    typedef u_int64_t Value;

    /** Constructor. */
    MPHF ()  {}

    /** Constructor. */
    MPHF (tools::dp::Iterator<Key>* iterator, size_t nbElts)  {  error();  }

    /** Constructor. */
    MPHF (tools::collections::Iterable<Key>* iterable)  { error(); }

    /** Build the hash function from a set of items.
     * \param[in] iterator : keys iterator
     * \param[in] nbItems : number of keys (if known) */
    void build (tools::dp::Iterator<Key>* iterator, size_t nbElts)   { error(); }

    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Value operator () (const Key& key) { error(); return Value(); }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { error(); return 0; }

private:

    /** Default error management. */
    void error () const { throw gatb::core::system::ExceptionNotImplemented(); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_MPHF_HPP_ */
