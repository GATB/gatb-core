/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Range.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterable interface
 */

#ifndef _GATB_CORE_TOOLS_MISC_RANGE_HPP_
#define _GATB_CORE_TOOLS_MISC_RANGE_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Tools package */
namespace tools     {
/** \brief Misc interfaces */
namespace misc      {
/********************************************************************************/

/** \brief Definition of an interval (inspired by std::pair). It is possible to define
 * 'reversed' range, ie. with a beginning greater than the end.
 */
template <class T> class Range
{
public:

    /** Default constructor. */
    Range() : begin(T()), end(T()) {}

    /** Constructor taking the [begin,end] couple as arguments.
     * \param[in] x : beginning of the range.
     * \param[in] y : end of the range.
     */
    Range(const T& x, const T& y) : begin(x), end(y) {}

    /** Copy constructor. */
    template <class U>  Range (const Range<U> &p) : begin(p.begin), end(p.end) { }

    /** \return the begin bound of the range. */
    T getBegin() { return begin; }

    /** \return the end bound of the range. */
    T getEnd  () { return end; }

    /** Returns the length of the range. */
    T getLength ()  const  { return (end >= begin ? end - begin + 1 : begin - end + 1); }

    /** Tells whether the provided value is included into the 'this' instance.
     * \param[in] val : value to be tested
     * \return true if the provided value is inside the current range, false otherwise.
     */
    bool includes (const T& val) const   {   return (this->begin <= val)  &&  (val <= this->end);  }

    /** Equality operator.
     * \param[in] r : the range to be compared to.
     * \return true if the ranges are the same (same beginning, same end), false otherwise. */
    bool operator== (const Range& r) const
    {
        return begin==r.begin && end==r.end;
    }

    /** InEquality operator.
     * \param[in] r : the range to be compared to.
     * \return false if the ranges are the same (same beginning, same end), true otherwise. */
    bool operator!= (const Range& r) const
    {
        return begin!=r.begin || end!=r.end;
    }

    /* */
    class Iterator : public dp::Iterator<T>
    {
    public:

        Iterator (const T& x, const T& y) : _begin(x), _end(y) {}

        Iterator (Range& ref) : _begin(ref.getBegin()), _end(ref.getEnd())  {}

        void first()  { *this->_item = _begin; }

        void next()   { (*this->_item) ++; }

        bool isDone() { return (*this->_item) > _end; }

        T& item ()     { return (*this->_item); }

    private:
        T      _begin;
        T      _end;
    };

    friend class Iterator;

private:
    /** The template (integer) type. */
    typedef T type;

    /** Beginning of the range. */
    T begin;

    /** End of the range. */
    T end;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_MISC_RANGE_HPP_ */
