#include <gatb/kmer/impl/PartitionsCommand.cpp>

/********************************************************************************/
namespace gatb { namespace core { namespace kmer { namespace impl  {
/********************************************************************************/
template class PartitionsCommand <KSIZE_5>;
template class PartitionsCommand <KSIZE_6>;
template class PartitionsCommand <KSIZE_7>;
template class PartitionsCommand <KSIZE_8>;

template class PartitionsByHashCommand <KSIZE_5>;
template class PartitionsByHashCommand <KSIZE_6>;
template class PartitionsByHashCommand <KSIZE_7>;
template class PartitionsByHashCommand <KSIZE_8>;

template class PartitionsByVectorCommand <KSIZE_5>;
template class PartitionsByVectorCommand <KSIZE_6>;
template class PartitionsByVectorCommand <KSIZE_7>;
template class PartitionsByVectorCommand <KSIZE_8>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
