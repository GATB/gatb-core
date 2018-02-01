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

#include <iostream>
#include <fstream>
#include <sstream>
#include <atomic>
#include <set>
#include <vector>
#include <string>
#include <mutex>
#include <gatb/tools/storage/impl/Storage.hpp>

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
            if (!_insertHandle) { std::cout << "error opening " << filename << " for writing." << std::endl; exit(1);}
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
            unsigned int res = fprintf (_insertHandle, "%s", buffer.c_str());
            if (res != buffer.size())            {  std::cout << "couldn't flush, written " << res << " out of " << buffer.size() << std::endl; exit(1);}
            if (_insertHandle    != 0)  { fflush  (_insertHandle);              }
            buffer_length = 0;
            buffer.clear();
        }
};

// not using BankFasta because I suspect that it does some funky memory fragmentation. so this one is unbuffered
class UnbufferedFastaIterator 
{
        std::ifstream *input;
    public:
        UnbufferedFastaIterator(const std::string &filename)
        {
            input = new std::ifstream(filename);
        }

        ~UnbufferedFastaIterator() { delete input;}

        bool read(std::string &seq, std::string &comment)
        {
            std::string line;
            if (std::getline(*input, line))
            {
                if (line.empty())
                    return false;
                if (line[0] == '>')
                    comment = line.substr(1);
            }
            else 
                return false;
            if (std::getline(*input, line))
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
            input->clear();
            input->seekg(0);
        }

};

    template<size_t SPAN>
void bglue(gatb::core::tools::storage::impl::Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int nb_glue_partitions, 
        int nb_threads, 
        bool verbose
        );

}}}}

#endif

