#include <gatb/kmer/impl/PartitionsCommand.cpp>

/********************************************************************************/
namespace gatb { namespace core { namespace kmer { namespace impl  {
/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class PartitionsCommand <KSIZE_1>;
template class PartitionsCommand <KSIZE_2>;
template class PartitionsCommand <KSIZE_3>;
template class PartitionsCommand <KSIZE_4>;

template class PartitionsByHashCommand <KSIZE_1>;
template class PartitionsByHashCommand <KSIZE_2>;
template class PartitionsByHashCommand <KSIZE_3>;
template class PartitionsByHashCommand <KSIZE_4>;

template class PartitionsByVectorCommand <KSIZE_1>;
template class PartitionsByVectorCommand <KSIZE_2>;
template class PartitionsByVectorCommand <KSIZE_3>;
template class PartitionsByVectorCommand <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
