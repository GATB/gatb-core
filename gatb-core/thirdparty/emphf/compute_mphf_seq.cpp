#include "compute_mphf_generic.hpp"
#include "hypergraph.hpp"
#include "hypergraph_sorter_seq.hpp"

int main(int argc, char** argv)
{
    using namespace emphf;
    return compute_mphf_main<hypergraph_sorter_seq<hypergraph<uint32_t>>,
                             hypergraph_sorter_seq<hypergraph<uint64_t>>,
                             jenkins64_hasher>(argc, argv);
}
