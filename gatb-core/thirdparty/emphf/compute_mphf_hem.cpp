#include <iostream>
#include <fstream>
#include <iterator>
#include <random>

#include "common.hpp"
#include "mphf_hem.hpp"
#include "mmap_memory_model.hpp"
#include "base_hash.hpp"
#include "perfutils.hpp"

int main(int argc, char** argv)
{
    using namespace emphf;

    if (argc < 2) {
        std::cerr << "Expected: " << argv[0]
                  << " <filename> [output_filename]" << std::endl;
        std::terminate();
    }

    const char* filename = argv[1];
    std::string output_filename;

    if (argc >= 3) {
        output_filename = argv[2];
    }

    logger() << "Processing " << filename << std::endl;

    file_lines lines(filename);
    size_t n = lines.size();
    logger() << n << " strings to process." << std::endl;

    stl_string_adaptor adaptor;
    mmap_memory_model mm;
    mphf_hem<jenkins32_hasher> mphf(mm, n, lines, adaptor);

    if (output_filename.size()) {
        std::ofstream os(output_filename, std::ios::binary);
        mphf.save(os);
    }

    return 0;
}
