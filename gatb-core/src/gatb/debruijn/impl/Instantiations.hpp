// i'm trying to reduce this file as much as I can, because nested templates are now nasty, and haven't found a way to properly specialize with NodeFast.


template<>   template <>  inline GraphT::Vector<NodeT> GraphT::successors   ( NodeT& node) const                 {  return getNodes(node, DIR_OUTCOMING); }
template<>   template <>  inline GraphT::Vector<NodeT> GraphT::predecessors ( NodeT& node) const                 {  return getNodes(node, DIR_INCOMING);  }

template<>   template <>  inline NodeT GraphT::neighbor ( NodeT& source, Direction dir, kmer::Nucleotide nt) const
{  bool exists=true; return getNode (source, dir, nt, exists);  }

template<>   template <>  inline NodeT GraphT::neighbor ( NodeT& source, Direction dir, kmer::Nucleotide nt, bool& exists) const
{  return getNode (source, dir, nt, exists);  }

template<>   template <>  inline NodeT GraphT::successor ( NodeT& source, kmer::Nucleotide nt) const
{  bool exists=true; return getNode (source, DIR_OUTCOMING, nt, exists);  }

template<>   template <>  inline NodeT GraphT::successor ( NodeT& source, kmer::Nucleotide nt, bool& exists) const
{  return getNode (source, DIR_OUTCOMING, nt, exists);  }

template<>   template <>  inline NodeT GraphT::predecessor ( NodeT& source, kmer::Nucleotide nt) const
{  bool exists=true; return getNode (source, DIR_INCOMING, nt, exists);  }

template<>   template <>  inline NodeT GraphT::predecessor ( NodeT& source, kmer::Nucleotide nt, bool& exists) const
{  return getNode (source, DIR_INCOMING, nt, exists);  }

template<>   template <>  inline GraphT::Vector<EdgeT> GraphT::successors   ( NodeT& node) const                 {  return getEdges(node, DIR_OUTCOMING); }
template<>   template <>  inline GraphT::Vector<EdgeT> GraphT::predecessors ( NodeT& node) const                 {  return getEdges(node, DIR_INCOMING);  }

/** so it appears that.. these strange two-nodes operations are used in Kissnp2 */
template<>   template <>  inline GraphT::Vector<std::pair<EdgeT,EdgeT> > GraphT::successors (const NodeT& node1, const NodeT& node2) const
{ return getEdgesCouple (node1, node2, DIR_OUTCOMING); }

template<>   template <>  inline GraphT::Vector<std::pair<NodeT,NodeT> > GraphT::successors (const NodeT& node1, const NodeT& node2) const
{ return getNodesCouple (node1, node2, DIR_OUTCOMING); }

template<>   template <>  inline GraphT::Vector<std::pair<EdgeT,EdgeT> > GraphT::predecessors (const NodeT& node1, const NodeT& node2) const
{ return getEdgesCouple (node1, node2, DIR_INCOMING); }

template<>   template <>  inline GraphT::Vector<std::pair<NodeT,NodeT> > GraphT::predecessors (const NodeT& node1, const NodeT& node2) const
{ return getNodesCouple (node1, node2, DIR_INCOMING); }

/** */
template<>   template <>  inline GraphT::Vector<BranchingNodeT> GraphT::successors (NodeT& node) const
{ return getBranchingNodeNeighbors (node, DIR_OUTCOMING);  }

template<>   template <>  inline GraphT::Vector<BranchingNodeT> GraphT::predecessors (NodeT& node) const
{ return getBranchingNodeNeighbors (node, DIR_INCOMING);  }

/** */
/** */
//template<>   template <> std::set<BranchingNodeT> GraphT::neighbors (std::set<BranchingNodeT>::iterator first, std::set<BranchingNode>::iterator last) const;
// later, cleanup if not used

/** */
template<>   template <>  inline GraphT::Vector<BranchingEdgeT> GraphT::successors (NodeT& node) const
{ return getBranchingEdgeNeighbors (node, DIR_OUTCOMING);  }

template<>   template <>  inline GraphT::Vector<BranchingEdgeT> GraphT::predecessors (NodeT& node) const
{ return getBranchingEdgeNeighbors (node, DIR_INCOMING);  }

/** */
template<>   template<> inline GraphT::Iterator<NodeT> GraphT::simplePath (NodeT& node, Direction dir) const  { return getSimpleNodeIterator(node, dir); }
template<>   template<> inline GraphT::Iterator<EdgeT> GraphT::simplePath (NodeT& node, Direction dir) const  { return getSimpleEdgeIterator(node, dir); }


