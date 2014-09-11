#include "compute_mphf_generic.hpp"
#include "mmap_memory_model.hpp"
#include "hypergraph_sorter_scan.hpp"

int main(int argc, char** argv)
{
    using namespace emphf;
    return compute_mphf_main<hypergraph_sorter_scan<uint32_t, mmap_memory_model>,
                             hypergraph_sorter_scan<uint64_t, mmap_memory_model>,
                             jenkins64_hasher>(argc, argv);
}
