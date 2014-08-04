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
class Histogram : public IHistogram, public system::SmartPointer
{
public:

    /** */
    Histogram (size_t length, tools::collections::Bag<Entry>* bag)
        : _length(length), _cutoff(0), _nbsolids(0),
          _histogram(0), _histogram_smoothed(0), _bag(0)
    {
        setBag (bag);

        _histogram = (Entry*) system::impl::System::memory().calloc (_length + 1, sizeof (Entry));
        memset (_histogram, 0, sizeof(Entry)*(_length + 1));

		_histogram_smoothed = (Entry*) system::impl::System::memory().calloc (_length + 1, sizeof (Entry));
        memset (_histogram_smoothed, 0, sizeof(Entry)*(_length + 1));
		
        for (size_t i=0; i<_length+1; i++)
        {
            _histogram[i].index     = i;
			_histogram_smoothed[i].index     = i;

            _histogram[i].abundance = 0;
        }
    }

    /** */
    virtual ~Histogram ()
    {
        setBag(0);
        system::impl::System::memory().free (_histogram);
        system::impl::System::memory().free (_histogram_smoothed);
    }

    /** */
    void inc (u_int16_t index)  { _histogram [(index >= _length) ? _length : index].abundance ++; }

    /** */
    void save ();

	
	void compute_threshold () ;
	u_int16_t get_solid_cutoff () ;
	u_int64_t get_nbsolids_auto () ;

    /** */
    size_t getLength() { return _length; }

    /** */
    u_int64_t& get (u_int16_t i)  { return _histogram[i].abundance; }

private:

    size_t  _length;
	u_int16_t _cutoff;
	u_int64_t _nbsolids;
	
    Entry*  _histogram;

	Entry*  _histogram_smoothed;

    tools::collections::Bag<Entry>* _bag;
    void setBag (tools::collections::Bag<Entry>* bag)  { SP_SETATTR(bag); }
};

/********************************************************************************/

/** */
class HistogramNull : public IHistogram, public system::SmartPointer
{
public:

    /** */
    void inc (u_int16_t abundance) {}

    /** */
    void save ()  {}
	
	u_int16_t get_solid_cutoff () {return 0; }
	u_int64_t get_nbsolids_auto () {return 0;}

	void compute_threshold () { }

    /** */
    size_t getLength() { return 0; }

    /** */
    u_int64_t& get (u_int16_t i)  { static u_int64_t foo; return foo; }
};

/********************************************************************************/

/** */
class HistogramCache : public IHistogram, public system::SmartPointer
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

	void compute_threshold () { return _ref->compute_threshold(); }

	u_int16_t get_solid_cutoff () {return _ref->get_solid_cutoff();}
	
	u_int64_t get_nbsolids_auto () {return _ref->get_nbsolids_auto();}

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
