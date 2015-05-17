#include <gatb/kmer/impl/DebloomMinimizerAlgorithm.cpp>

/********************************************************************************/
namespace gatb { namespace core { namespace kmer { namespace impl  {
/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomMinimizerAlgorithm <KSIZE_5>;
template class DebloomMinimizerAlgorithm <KSIZE_6>;
template class DebloomMinimizerAlgorithm <KSIZE_7>;
template class DebloomMinimizerAlgorithm <KSIZE_8>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
