/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file BankStrings.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Hard coded genomic bank (mainly for tests)
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_STRINGS_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_STRINGS_HPP_

/********************************************************************************/

#include <gatb/bank/impl/AbstractBank.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <vector>
#include <string>
#include <sstream>
#include <stdarg.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Interface for reading genomic databases.
 */
class BankStrings : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "strings"; }

    /** Constructor. */
    BankStrings (const std::vector<std::string>& sequencesData)
        : _sequencesData(sequencesData), _totalSize(0), _maxSize(0)  {  init ();  }

    /** Constructor. */
    BankStrings (const char* sequencesData[], size_t nb)  : _totalSize(0), _maxSize(0)
    {
        for (size_t i=0; i<nb; i++)  { _sequencesData.push_back (sequencesData[i]); }
        init ();
    }

    /** Constructor. */
    BankStrings (const char* seq, ...) : _totalSize(0), _maxSize(0)
    {
        va_list ap;
        va_start (ap, seq);
        for (const char* loop=seq; loop != 0; loop = va_arg (ap, const char*))  { _sequencesData.push_back(loop); }
        va_end   (ap);

        init ();
    }

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return new tools::dp::impl::VectorIterator2<Sequence> (_sequences); }

    /** \copydoc IBank::getNbItems */
    int64_t getNbItems () { return _sequences.size(); }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item)
    {
        _totalSize += item.getDataSize();
        if (_maxSize < item.getDataSize())  {  _maxSize = item.getDataSize(); }
        _sequences.push_back (item);
    }

    /** \copydoc IBank::flush */
    void flush ()  {}

    /** \copydoc IBank::getSize */
    u_int64_t getSize ()  { return _totalSize; }

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
    {
        number    = _sequences.size();
        totalSize = _totalSize;
        maxSize   = _maxSize;
    }

protected:

    std::vector<std::string> _sequencesData;

    /** Sequences. */
    std::vector<Sequence> _sequences;

    u_int64_t _totalSize;
    u_int64_t _maxSize;

    /** */
    void init ()
    {
        _sequences.resize (_sequencesData.size());

        for (size_t i=0; i<_sequences.size(); i++)
        {
            /** Shortcut. */
            Sequence& s = _sequences[i];

            s.setIndex (i);
            std::stringstream ss;  ss << "seq_" << i;   s._comment = ss.str();
            s.getData().setRef ((char*) _sequencesData[i].data(), _sequencesData[i].size());

            _totalSize += _sequencesData[i].size();

            if (_maxSize < _sequencesData[i].size())  {  _maxSize = _sequencesData[i].size(); }
        }
    }

};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_STRINGS_HPP_ */
