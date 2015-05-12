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

#ifndef _COUNT_PROCESSOR_CHAIN_HPP_
#define _COUNT_PROCESSOR_CHAIN_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <cstdarg>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorChain : public CountProcessorAbstract<span>
{
public:

    typedef ICountProcessor<span> CountProcessor;
    typedef typename Kmer<span>::Type           Type;

    CountProcessorChain (const std::vector<CountProcessor*> items) : _items(items) {}

    CountProcessorChain (CountProcessor* first, ...)
    {
        va_list ap;
        va_start (ap, first);
        for (CountProcessor* loop = first; loop != 0; loop = va_arg(ap, CountProcessor*))
        {
            _items.push_back (loop);  loop->use();
        }
        va_end (ap);
    }

    virtual ~CountProcessorChain()  {  for (size_t i=0; i<_items.size(); i++)  {  _items[i]->forget(); }  }

    /** */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0)
    {
        if (sum==0)  { sum = this->computeSum(count); }
        bool res = true;
        for (size_t i=0; res && i<_items.size(); i++)  {  res = _items[i]->process (partId, kmer, count, sum);  }
        return res;
    }

    /** */
    void begin (const Configuration& config)  {  for (size_t i=0; i<_items.size(); i++)  {  _items[i]->begin (config);  }  }
    void end   ()  {  for (size_t i=0; i<_items.size(); i++)  {  _items[i]->end   ();         }  }

    /** */
    void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name)
    {
        for (size_t i=0; i<_items.size(); i++)  {  _items[i]->beginPart (passId, partId, cacheSize, name);  }
    }

    void endPart   (size_t passId, size_t partId)
    {
        for (size_t i=0; i<_items.size(); i++)  {  _items[i]->endPart   (passId, partId);  }
    }

    /** */
    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;
        for (size_t i=0; i<_items.size(); i++)  {  result.add (0, _items[i]->getProperties());  }
        return result;
    }

    /** */
    std::vector<CountProcessor*> getInstances () const
    {
        std::vector<CountProcessor*> res;
        for (size_t i=0; i<_items.size(); i++)
        {
            std::vector<CountProcessor*> c = _items[i]->getInstances();
            for (size_t i=0; i<c.size(); i++)  { res.push_back(c[i]); }
        }
        return res;
    }


protected:

    /** */
    CountProcessorAbstract<span>* doClone ()
    {
        std::vector<CountProcessor*> clones;
        for (size_t i=0; i<_items.size(); i++)  {  clones.push_back (_items[i]->clone());  }
        return new CountProcessorChain (clones);
    }

    std::vector<CountProcessor*> _items;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_CHAIN_HPP_ */
