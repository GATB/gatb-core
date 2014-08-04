#include <iostream>
#include <fstream>

#include "test_mphf_generic.hpp"
#include "mphf.hpp"
#include "base_hash.hpp"

int main(int argc, char** argv)
{
    using namespace emphf;
    return test_mphf_main<mphf<jenkins64_hasher>>(argc, argv);
}
