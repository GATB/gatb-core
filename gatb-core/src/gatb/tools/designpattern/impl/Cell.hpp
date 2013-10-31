/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Node.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of INode interface
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_CELL_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_CELL_HPP_

#include <gatb/tools/designpattern/api/ICell.hpp>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** \brief Partial implementation of the INode interface
 *
 * The 'remove' method is still abstract.
 */
class Cell : public virtual ICell, public system::SmartPointer
{
public:

    /** Constructor. */
    Cell (ICell* parent, const std::string& id)  : _parent(0), _id(id)  { setParent(parent); }

    /** Destructor. */
    ~Cell ()   {  setParent(0);  }

    /** \copydoc ICell::getParent  */
    ICell* getParent () const { return _parent; }

    /** \copydoc ICell::getId  */
    const std::string& getId ()  const { return _id; }

    /** \copydoc ICell::getFullId  */
    std::string getFullId (char sep='.') const
    {
        if (_parent != 0)   {  return _parent->getId() + sep + getId();  }
        else                {  return getId();  }
    }

private:

    ICell* _parent;
    void setParent (ICell* parent)
    {
        _parent = parent;
        //SP_SETATTR(parent);
    }

    std::string _id;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_CELL_HPP_ */
