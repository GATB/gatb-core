/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/** \file Node.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Node
 */

#ifndef _GATB_CORE_DP_ICELL_HPP_
#define _GATB_CORE_DP_ICELL_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/SmartPointer.hpp>
#include <string>

/********************************************************************************/
namespace gatb  {
namespace core  {
/** \brief Tools package */
namespace tools {
/** \brief Design Patterns concepts */
namespace dp    {
/********************************************************************************/

class ICell : public virtual dp::ISmartPointer
{
public:

    /** Destructor. */
    virtual ~ICell() {}

    /** Get the parent node (if any).
     * \return the parent node. */
    virtual ICell* getParent () const = 0;

    /** Return the identifier of the node.
     * \return the id of the node. */
    virtual const std::string& getId () const = 0;

    /** Return the full identifier (like a path "x.y.z")
     * \return the full identifier. */
    virtual std::string getFullId () const = 0;

    /** Physically remove the node. */
    virtual void remove () = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ICELL_HPP_ */
