#include <gatb/kmer/impl/SortingCountAlgorithm.cpp>

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html
// (last example)
// also, to reduce compilation time, I'm splitting it into 4 files that will be compiled in parallel

/********************************************************************************/
namespace gatb { namespace core { namespace kmer { namespace impl  {
/********************************************************************************/

template class SortingCountAlgorithm <KSIZE_1>;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
