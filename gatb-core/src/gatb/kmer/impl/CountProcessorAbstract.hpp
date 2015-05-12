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

#ifndef _COUNT_PROCESSOR_ABSTRACT_HPP_
#define _COUNT_PROCESSOR_ABSTRACT_HPP_

/********************************************************************************/

#include <gatb/kmer/api/ICountProcessor.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** Abstract implementation of ICountProcessor interface.
 */
template<size_t span>
class CountProcessorAbstract : public ICountProcessor<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type Type;

    /** Constructor. */
    CountProcessorAbstract ()  : _prototype(0) {}

    /** Destructor. */
    virtual ~CountProcessorAbstract ()  {}

    /** \copydoc ICountProcessor<span>::clone */
    ICountProcessor<span>* clone ()
    {
        CountProcessorAbstract<span>* result = this->doClone();
        result->setPrototype (this);
        return result;
    }

    /** \copydoc ICountProcessor<span>::process */
    virtual bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0)  {  return true;  }

    /** \copydoc ICountProcessor<span>::begin */
    virtual void begin(const Configuration& config) {}

    /** \copydoc ICountProcessor<span>::end */
    virtual void end  () {}

    /** \copydoc ICountProcessor<span>::beginPart */
    virtual void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name) {}

    /** \copydoc ICountProcessor<span>::endPart */
    virtual void endPart   (size_t passId, size_t partId) {}

    /** \copydoc ICountProcessor<span>::getProperties */
    virtual tools::misc::impl::Properties getProperties() const { return tools::misc::impl::Properties(); }

    /** \copydoc ICountProcessor<span>::getInstances */
    virtual std::vector<ICountProcessor<span>*> getInstances () const
    {
        std::vector<ICountProcessor<span>*> res;
        res.push_back((CountProcessorAbstract*)this);
        return res;
    }

protected:

    /** Abstract method to be implemented by children class. */
    virtual CountProcessorAbstract* doClone () = 0;

    /** \copydoc ICountProcessor<span>::computeSum */
    CountNumber computeSum (const CountVector& count) const
    {
        /** Optimization. */
        if (count.size()==1)  { return count[0]; }
        CountNumber sum=0; for (size_t k=0; k<count.size(); k++)  { sum+=count[k]; }  return sum;
    }

    /** Getter on the prototype (if any)
     * \return the prototype instance. */
    ICountProcessor<span>* getPrototype ()  { return _prototype; }

private:

    ICountProcessor<span>* _prototype;
    void setPrototype (ICountProcessor<span>* proto)  { _prototype = proto; }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_ABSTRACT_HPP_ */
