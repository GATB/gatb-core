/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Histogram.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Histogram feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_
#define _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/misc/api/IHistogram.hpp>
#include <string>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** TBD */
class Histogram : public IHistogram, public dp::SmartPointer
{
public:

    /** */
    Histogram (size_t length, tools::collections::Bag<Entry>* bag) : _length(length), _bag(0)
    {
        setBag (bag);

        _histogram = (Entry*) system::impl::System::memory().calloc (_length + 1, sizeof (Entry));
        memset (_histogram, 0, sizeof(Entry)*(_length + 1));

        for (size_t i=0; i<_length+1; i++)
        {
            _histogram[i].index     = i;
            _histogram[i].abundance = 0;
        }
    }

    /** */
    virtual ~Histogram ()
    {
        setBag(0);
        system::impl::System::memory().free (_histogram);
    }

    /** */
    void inc (u_int16_t index)  { _histogram [(index >= _length) ? _length : index].abundance ++; }

    /** */
    void save ();

    /** */
    size_t getLength() { return _length; }

    /** */
    u_int64_t& get (u_int16_t i)  { return _histogram[i].abundance; }

private:

    size_t  _length;
    Entry*  _histogram;

    tools::collections::Bag<Entry>* _bag;
    void setBag (tools::collections::Bag<Entry>* bag)  { SP_SETATTR(bag); }
};

/********************************************************************************/

/** */
class HistogramNull : public IHistogram, public dp::SmartPointer
{
public:

    /** */
    void inc (u_int32_t abundance) {}

    /** */
    void save ()  {}

    /** */
    size_t getLength() { return 0; }

    /** */
    u_int64_t& get (u_int16_t i)  { static u_int64_t foo; return foo; }
};

/********************************************************************************/

/** */
class HistogramCache : public IHistogram, public dp::SmartPointer
{
public:

    /** */
    HistogramCache (IHistogram* ref, system::ISynchronizer* synchro=0)
        : _ref(0), _synchro(synchro), _localHisto(ref ? ref->getLength() : 0, 0) {  setRef(ref); }

    /** */
    ~HistogramCache()
    {
        system::LocalSynchronizer ls (_synchro);
        for (size_t cc=1; cc<_localHisto.getLength(); cc++)  {  _ref->get(cc) += _localHisto.get(cc);  }
        setRef (0);
    }

    /** */
    void inc (u_int16_t index)  { _localHisto.inc (index); }

    /** */
    void save ()  { return _ref->save(); }

    /** */
    size_t getLength() { return _localHisto.getLength(); }

    /** */
    u_int64_t& get (u_int16_t i)  { return _localHisto.get(i); }

private:

    IHistogram* _ref;
    void setRef (IHistogram* ref)  { SP_SETATTR(ref); }

    system::ISynchronizer* _synchro;
    Histogram              _localHisto;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_ */
