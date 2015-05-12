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

#ifndef _COUNT_PROCESSOR_HISTOGRAM_HPP_
#define _COUNT_PROCESSOR_HISTOGRAM_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <cstdarg>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorHistogram : public CountProcessorAbstract<span>
{
public:

    typedef typename Kmer<span>::Type Type;

    CountProcessorHistogram (
        tools::storage::impl::Group& group,
        size_t histoMax           = 10000,
        size_t min_auto_threshold = 3
    )
        : _group(group), _histogram(0), _min_auto_threshold(min_auto_threshold)
    {
        setHistogram (new tools::misc::impl::Histogram (histoMax));
    }

    CountProcessorHistogram (
        tools::storage::impl::Group& group,
        tools::misc::IHistogram* histogram,
        size_t min_auto_threshold = 3
    )
        : _group(group), _histogram(0), _min_auto_threshold(min_auto_threshold)
    {
        setHistogram (histogram);
    }

    virtual ~CountProcessorHistogram ()
    {
        setHistogram(0);
    }

    /** */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
        _histogram->inc (sum);
        return true;
    }

    void end ()
    {
        using namespace tools::math;

        /** We save the histogram if any. */
        _histogram->save (_group);

        /** compute auto cutoff **/
        _histogram->compute_threshold (_min_auto_threshold);

        /** store auto cutoff and corresponding number of solid kmers **/
        tools::collections::Collection<NativeInt64>& storecutoff = _group.getCollection<NativeInt64>("cutoff") ;
        storecutoff.insert(_histogram->get_solid_cutoff());
        storecutoff.flush();

        tools::collections::Collection<NativeInt64>& storesolids = _group.getCollection<NativeInt64>("nbsolidsforcutoff") ;
        storesolids.insert(_histogram->get_nbsolids_auto());
        storesolids.flush();
    }

    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;

        result.add (0, "histogram");
        result.add (1, "cutoff",            "%ld",  _histogram->get_solid_cutoff());
        result.add (1, "nb_ge_cutoff",      "%ld",  _histogram->get_nbsolids_auto());
        // result->add (1, "percent_ge_cutoff", "%.1f", nbSolids > 0 ? 100.0 * (double)_histogram->get_nbsolids_auto() / (double)_bankStats.kmersNbValid : 0);
        result.add (1, "first_peak",         "%ld",  _histogram->get_first_peak());

        // double N = ((double)_histogram->get_first_peak() * _bankStats.getSeqMean()) / (_bankStats.getSeqMean() - _kmerSize + 1);
        // if (N > 0)  {  getInfo()->add (3, "genome_size_estimate", "%.0f",  (double)_bankStats.sequencesTotalLength / N);  }

        return result;
    }

    /** */
    gatb::core::tools::misc::IHistogram* getHistogram() { return _histogram; }


protected:

    /** */
    CountProcessorAbstract<span>* doClone ()
    {
        /** We encapsulate the histogram with a cache. */
        return new CountProcessorHistogram (
            _group,
            new gatb::core::tools::misc::impl::HistogramCache (_histogram),
            _min_auto_threshold
        );
    }

private:

    tools::storage::impl::Group& _group;

    gatb::core::tools::misc::IHistogram* _histogram;
    void setHistogram (gatb::core::tools::misc::IHistogram* histogram)  { SP_SETATTR(histogram); }

    size_t _min_auto_threshold;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_HISTOGRAM_HPP_ */
