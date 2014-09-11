#include <iostream>
#include <fstream>

#include "test_mphf_generic.hpp"
#include "mphf_hem.hpp"
#include "base_hash.hpp"

int main(int argc, char** argv)
{
    using namespace emphf;
    return test_mphf_main<mphf_hem<jenkins32_hasher>>(argc, argv);
}
