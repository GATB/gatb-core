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

/** \file ContainerNode.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Container implementation
 */

#ifndef _GATB_CORE_DEBRUIJN_IMPL_CONTAINER_NODE_HPP_
#define _GATB_CORE_DEBRUIJN_IMPL_CONTAINER_NODE_HPP_

/********************************************************************************/

#include <gatb/debruijn/api/IContainerNode.hpp>
#include <cstdarg>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace debruijn  {
namespace impl      {
/********************************************************************************/

/** \brief IContainerNode implementation with a Bloom filter and a cFP set
 *
 *  In the GATB terminology, this object contains the information relative to the nodes of the dBG.
 *  It is not a set of nodes, as we don't store nodes explicitly. 
 *  Only one operation is supported:
 *     contains()
 *
 *  In the ContainerNode implementation, this object is actually just the Bloom filter + the set of False positives.
 */
template <typename Item> class ContainerNode : public IContainerNode<Item>, public system::SmartPointer
{
public:

    /** Constructor
     * \param[in] bloom : the Bloom filter.
     * \param[in] falsePositives : the cFP container. */
    ContainerNode (
        tools::collections::Container<Item>* bloom,
        tools::collections::Container<Item>* falsePositives
	) : _bloom(0), _falsePositives(0)
    {
		setBloom(bloom);
		setFalsePositives(falsePositives);
    }

	/** Destructor. */
	~ContainerNode ()
	{
		setBloom(0);
		setFalsePositives(0);
	}

    /** \copydoc IContainerNode::contains */
    bool contains (const Item& item)  {  return (_bloom->contains(item) && !_falsePositives->contains(item));  }

protected:

    tools::collections::Container<Item>* _bloom;
	void setBloom (tools::collections::Container<Item>* bloom)  { SP_SETATTR(bloom); }

	tools::collections::Container<Item>* _falsePositives;
	void setFalsePositives (tools::collections::Container<Item>* falsePositives)  { SP_SETATTR(falsePositives); }
};

/********************************************************************************/

/** \brief IContainerNode implementation with a Bloom filter
 *
 * This implementation has no critical False Positive set, so it implies that
 * Graph instances using it will have false positive nodes (old 'Titus' mode).
 */
template <typename Item> class ContainerNodeNoCFP : public ContainerNode<Item>
{
public:

    /** Constructor
     * \param[in] bloom : the Bloom filter. */
    ContainerNodeNoCFP (tools::collections::Container<Item>* bloom) : ContainerNode<Item>(bloom, NULL) {}

    /** \copydoc IContainerNode::contains */
    bool contains (const Item& item)  {  return (this->_bloom)->contains(item);  }
};

/********************************************************************************/
/** \brief IContainerNode implementation with cascading Bloom filters
 *
 * This implementation uses cascading Bloom filters for coding the cFP set.
 */
template <typename Item> class ContainerNodeCascading : public IContainerNode<Item>, public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] bloom : the Bloom filter.
     * \param[in] bloom2 : first Bloom filter of the cascading Bloom filters
     * \param[in] bloom3 : second Bloom filter of the cascading Bloom filters
     * \param[in] bloom4 : third Bloom filter of the cascading Bloom filters
     * \param[in] falsePositives : false positives container
     */
    ContainerNodeCascading (
        tools::collections::Container<Item>* bloom,
        tools::collections::Container<Item>* bloom2,
        tools::collections::Container<Item>* bloom3,
        tools::collections::Container<Item>* bloom4,
        tools::collections::Container<Item>* falsePositives
	) : _bloom(0), _bloom2(0), _bloom3(0), _bloom4(0), _falsePositives(0), _cfpArray(4)
    {
		setBloom          (bloom);
        setBloom2         (bloom2);
		setBloom3         (bloom3);
		setBloom4         (bloom4);
		setFalsePositives (falsePositives);

		_cfpArray[0] = bloom2;
        _cfpArray[1] = bloom3;
        _cfpArray[2] = bloom4;
        _cfpArray[3] = falsePositives;
    }

	/** Destructor */
	~ContainerNodeCascading ()
	{
        setBloom          (0);
		setBloom2         (0);
		setBloom3         (0);
		setBloom4         (0);
		setFalsePositives (0);
	}

    /** \copydoc IContainerNode::contains */
    bool contains (const Item& item)  {  return (_bloom->contains(item) && ! containsCFP(item));  }

private:

    tools::collections::Container<Item>* _bloom;
    void setBloom (tools::collections::Container<Item>* bloom)  { SP_SETATTR(bloom); }

    tools::collections::Container<Item>* _bloom2;
	void setBloom2 (tools::collections::Container<Item>* bloom2)  { SP_SETATTR(bloom2); }

	tools::collections::Container<Item>* _bloom3;
	void setBloom3 (tools::collections::Container<Item>* bloom3)  { SP_SETATTR(bloom3); }

	tools::collections::Container<Item>* _bloom4;
	void setBloom4 (tools::collections::Container<Item>* bloom4)  { SP_SETATTR(bloom4); }

	tools::collections::Container<Item>* _falsePositives;
	void setFalsePositives (tools::collections::Container<Item>* falsePositives)  { SP_SETATTR(falsePositives); }

    std::vector<tools::collections::Container<Item>*> _cfpArray;

    /** */
    bool containsCFP (const Item& item)
    {
        if (_bloom2->contains(item))
        {
            if (!_bloom3->contains(item))
                return true;

            else if (_bloom4->contains(item) && !_falsePositives->contains(item))
                return true;
        }
        return false;
    }

};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DEBRUIJN_IMPL_CONTAINER_NODE_HPP_ */
