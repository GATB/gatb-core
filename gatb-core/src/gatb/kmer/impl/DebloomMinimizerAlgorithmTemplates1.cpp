#include <gatb/kmer/impl/DebloomMinimizerAlgorithm.cpp>

/********************************************************************************/
namespace gatb { namespace core { namespace kmer { namespace impl  {
/********************************************************************************/

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

template class DebloomMinimizerAlgorithm <KSIZE_1>;
template class DebloomMinimizerAlgorithm <KSIZE_2>;
template class DebloomMinimizerAlgorithm <KSIZE_3>;
template class DebloomMinimizerAlgorithm <KSIZE_4>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
