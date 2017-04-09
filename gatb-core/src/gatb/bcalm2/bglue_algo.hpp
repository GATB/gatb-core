/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014-2016  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/


#ifndef _GATB_CORE_BGLUE_ALGO_HPP_
#define _GATB_CORE_BGLUE_ALGO_HPP_

#include "unionFind.hpp"
#include <atomic>
#include <set>
#include <vector>
#include <string>
#include <mutex>
#include <unordered_map>
#include <BooPHF/BooPHF.h>
#include <ctime> // for time
#include <iostream> // for time (and maybe other things?)
#include <iomanip> // for cout mods
#include "ThreadPool.h"
/*#include "ctpl_stl.h" // alternative to threadpool // https://github.com/vit-vit/CTPL/blob/master/ctpl_stl.h // didn't commit because didnt use
#include "buffer_allocator.h" // memory pool from https://github.com/vincetse/allocator, didn't commit the files because didnt use
#include "buffer_manager.h" // memory pool
*/

#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/storage/impl/StorageTools.hpp>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/kmer/impl/PartiInfo.hpp>   // for repartitor 
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/collections/impl/BooPHF.hpp>


//heh at this point I could have maybe just included gatb_core.hpp but well, no circular dependencies, this file is part of gatb-core now.

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;



namespace gatb { namespace core { namespace debruijn { namespace impl  {

// buffered and also thread-safe thank to a lock
// not using BankFasta because I dont want to be recording variable-length strings in a std::vector<>, potential memory fragmentation
// so instead it's in a flat buffer
class BufferedFasta
{
        std::mutex mtx;
        std::string buffer;
        unsigned long buffer_length;
        FILE* _insertHandle;

    public:
        unsigned long max_buffer;
        bool threadsafe;

        BufferedFasta(const std::string filename, unsigned long given_max_buffer = 50000)
        {
            max_buffer = given_max_buffer; // that much of buffering will be written to the file at once (in bytes)
            threadsafe = true;
            buffer_length = 0;
            _insertHandle = fopen (filename.c_str(), "w");
            buffer.reserve(max_buffer+1000/*security*/);
        }

        ~BufferedFasta()
        {
            flush();
            fclose(_insertHandle);
            std::string().swap(buffer);
        }

        void insert(const std::string &seq, const std::string &comment)
        {
            if (threadsafe)
                mtx.lock();
            unsigned int insert_size = seq.size() + comment.size() + 3;
            if (buffer_length + insert_size > max_buffer)
                flush();
            buffer_length += insert_size;
            buffer += ">" + comment + "\n" + seq + "\n";
            if (threadsafe)
                mtx.unlock();
        }

        void flush()
        {
            fprintf (_insertHandle, "%s", buffer.c_str());
            if (_insertHandle    != 0)  { fflush  (_insertHandle);              }
            buffer_length = 0;
            buffer.clear();
        }
};

// not using BankFasta because I suspect that it does some funky memory fragmentation. so this one is unbuffered
class UnbufferedFastaIterator 
{
        std::ifstream input;
    public:
        UnbufferedFastaIterator(const std::string &filename)
        {
            std::ifstream(filename).swap(input);
        }

        ~UnbufferedFastaIterator() {}

        bool read(std::string &seq, std::string &comment)
        {
            std::string line;
            if (std::getline(input, line))
            {
                if (line.empty())
                    return false;
                if (line[0] == '>')
                    comment = line.substr(1);
            }
            else 
                return false;
            if (std::getline(input, line))
            {
                if (line.empty())
                {
                    std::cout << "unexpected entry of one-line fasta: " <<  comment << std::endl;
                    exit(1);
                }
                seq = line;
                return true;
            }
            return false;
        }

        void restart()
        {
            input.clear();
            input.seekg(0);
        }

};

    template<size_t SPAN>
void bglue(Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int nb_threads, 
        bool verbose
        );

}}}}

#endif

