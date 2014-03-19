/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
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

#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <gatb/system/impl/System.hpp>

#include <errno.h>

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;

#define DEBUG(a)  //printf a

#define BINREADS_BUFFER 100000

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

int NT2int(char nt)
{
    int i;
    i = nt;
    i = (i>>1)&3; // that's quite clever, guillaume.
    return i;
}

unsigned char  code4NT(const char *seq)
{
    int i;
    unsigned char x;

    x=0;
    for (i=0; i<4; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}

unsigned char  code_n_NT(const char *seq, int nb)
{
    int i;
    unsigned char x;

    x=0;
    for (i=0; i<nb; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    x = x << ((4-nb)*2) ;
    return x;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankBinary::BankBinary (const std::string& filename, size_t nbValidLetters)
    : _filename(filename), _nbValidLetters(nbValidLetters), binary_read_file(0)
{
    read_write_buffer_size = BINREADS_BUFFER;

    //open (true);

    buffer = (unsigned char *) malloc(read_write_buffer_size*sizeof(unsigned char));

    cpt_buffer = 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankBinary::~BankBinary ()
{
    if(buffer!=NULL)
    {
        free (buffer); //buffer =NULL;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankBinary::insert (const Sequence& seq)
{
    /** Shortcuts. */
    int   whole_readlen = seq.getDataSize();
    char* pt_start      = seq.getDataBuffer();
    char* pt_end        = pt_start + whole_readlen;
    
    int readlen;
    int tai = readlen;
    unsigned char rbin;
    unsigned int block_size = 0;
    char *pt;
    
    char * pt_begin = pt_start;
    int idx =0 ;

    for ( ; pt_begin < pt_end; pt_begin += idx)
    {
        idx=0; // start a new read

        if (_nbValidLetters > 0)
        {
            //skips NN
            while (pt_begin < pt_end && *pt_begin =='N')
            {
                pt_begin ++;
            }
            // goes to next N or end of seq
            while ( ((pt_begin +idx) < pt_end) && (pt_begin[idx] !='N') )
            {
                idx++;
            }
        }
        else
        {
            idx = whole_readlen;
        }

        if (0<_nbValidLetters && idx<_nbValidLetters)  {  continue; }

        //we have a seq beginning at  pt_begin of size idx  ,without any N, will be treated as a read: of size readlen, beginning at pt
        readlen = tai = idx;
        pt      = pt_begin;
        
        /** We may have to open the file at first call. */
        if (binary_read_file == 0)  {  open (true); }
        
        //todo : also flush to disk  sometimes (ie if very large buffer, to create smaller blocks..)
        if(cpt_buffer >= (read_write_buffer_size-readlen) || cpt_buffer > 10000000 )  ////not enough space to store next read   true space is 4 + readlen/4 + rem
            //flush buffer to disk
        {
            block_size = cpt_buffer;
            
            if(block_size) fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
            if (!fwrite(buffer, 1, cpt_buffer, binary_read_file)) // write a block, it ends at end of a read
            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }
            cpt_buffer=0;
        }
        
        //check if still not enough space in empty buffer : can happen if large read, then enlarge buffer
        if(read_write_buffer_size < readlen)
        {
            read_write_buffer_size = 2*readlen; // too large but ok
            buffer =  (unsigned char *) realloc(buffer,sizeof(unsigned char) * read_write_buffer_size);
        }
        
        /** We write the length of the read. */
        memcpy(buffer+cpt_buffer,&readlen,sizeof(int));
        cpt_buffer+= sizeof(int);
        
        /** We write one byte for 4 nucleotides. */
        for (tai=readlen; tai>=4  ; tai-=4)
        {
            rbin = code4NT(pt);
            buffer[cpt_buffer]=rbin; cpt_buffer++;
            pt +=4;
        }
        
        /** We write the remaining part. */
        if(tai)
        {
            rbin = code_n_NT(pt,tai);
            buffer[cpt_buffer]=rbin; cpt_buffer++;
        }
    }
 }

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankBinary::flush ()
{
    unsigned int block_size =0;

    if (binary_read_file != 0)
    {
        if(cpt_buffer)
        {
            block_size = cpt_buffer;

            fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
            if (!fwrite(buffer, 1, cpt_buffer, binary_read_file))
            {
                throw gatb::core::system::ExceptionErrno (STR_BANK_unable_write_file);
            }
        }
        cpt_buffer=0;

        fclose(binary_read_file);
        binary_read_file = 0;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankBinary::open (bool write)
{
    binary_read_file = fopen (_filename.c_str(), write?"wb":"rb");
    if( binary_read_file == NULL)
    {
        throw gatb::core::system::ExceptionErrno (STR_BANK_unable_open_file, _filename.c_str());
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
u_int64_t BankBinary::getSize ()
{
    return System::file().getSize (_filename);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankBinary::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    /** We create an iterator for the bank. */
    BankBinary::Iterator it (*this);

    /** We return the estimation of sequences information. */
    it.estimate (number, totalSize, maxSize);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankBinary::Iterator::Iterator (BankBinary& ref)
    : _ref(ref), _isDone(true), cpt_buffer(0), blocksize_toread(0), nseq_lues(0),
      read_write_buffer_size(BINREADS_BUFFER), binary_read_file(0),
      _bufferData (0), _index(0)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankBinary::Iterator::~Iterator ()
{
    if (binary_read_file != 0)  {  fclose (binary_read_file);  }

    setBufferData(0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankBinary::Iterator::first()
{
    if (binary_read_file == 0)
    {
        /** We open the binary file at first call. */
        binary_read_file = fopen (_ref._filename.c_str(), "rb");

        /** We check that the file is opened => send an exception otherwise. */
        if (binary_read_file == 0)  {  throw gatb::core::system::ExceptionErrno (STR_BANK_unable_open_file, _ref._filename.c_str());  }
    }

    if (binary_read_file != 0)  {  rewind (binary_read_file); }

    /** We reinitialize some attributes. */
    _isDone          = false;
    cpt_buffer       = 0;
    blocksize_toread = 0;
    nseq_lues        = 0;

    /** We go to the next sequence. */
    next();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankBinary::Iterator::next ()
{
    int len = 0;
    unsigned int block_size = 0;

    //////////////////////////////////////////////
    //reading new block from disk if needed
    //////////////////////////////////////////////
    if (cpt_buffer == blocksize_toread)
    {
        /** We read the size of the following cache buffer. */
        if (! fread(&block_size,sizeof(unsigned int),1, binary_read_file)) //read block header
        {
            _isDone = true;
            return;
        }

        /** We are about to read another chunk of data from the disk. We need */
        setBufferData (new Data (block_size));

        size_t nbRead = fread (_bufferData->getBuffer(), sizeof( char),block_size, binary_read_file); // read a block of reads into the buffer

        cpt_buffer       = 0;
        blocksize_toread = block_size;
    }

    //////////////////////////////////////////////
    //now parse the whole block in ram
    //////////////////////////////////////////////
    if (cpt_buffer < blocksize_toread)
    {
        /** We increase the number of read sequences so far. */
        nseq_lues ++;

        memcpy (&len, _bufferData->getBuffer() + cpt_buffer, sizeof(int)); // read len

        /** We go ahead in the file parsing. */
        cpt_buffer += sizeof(int);

        int nchar = (len+3)/4;

        /** We update the information of the current sequence.
         * NOTE: we keep the original size of the data, not the compressed one. */
        _item->setDataRef (_bufferData, cpt_buffer, len);

        /** We set the sequence index. */
        _item->setIndex (_index++);

        /** We go ahead in the file parsing. */
        cpt_buffer += nchar;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void  BankBinary::Iterator::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    /** We initialize the provided arguments. */
    number    = 0;
    totalSize = 0;
    maxSize   = 0;

    /** We open the binary file at first call. */
    IFile* file = System::file().newFile (_ref._filename, "rb");
    if (file != 0)
    {
        vector<char> buffer;

        while (file->isEOF() == false)
        {
            unsigned int block_size = 0;

            // we read the block header
            if (! file->fread (&block_size,sizeof(unsigned int), 1))  { break; }

            // we resize the buffer
            buffer.resize (block_size);

            // read a block of reads into the buffer
            if (! file->fread (buffer.data(), sizeof(char), block_size))  { break; }

            // we loop the sequences
            char* loop = buffer.data();
            while (loop < buffer.data() + block_size)
            {
                int readlen = 0;
                memcpy (&readlen, loop, sizeof(int));

                loop += sizeof(int) + (readlen+3)/4;

                number ++;
                if (readlen > maxSize)  { maxSize = readlen; }
                totalSize += readlen;
            }

            // we may count only a part of the file
            if (number >= _ref.getEstimateThreshold())  { break; }
        }

        // we extrapolate the result according to the current location in the file
        if (file->isEOF() == false)
        {
            // we keep the current location in the file
            u_int64_t current = file->tell ();

            // we go to the end of the file
            file->seeko (0, SEEK_END);

            // we keep the current location in the file
            u_int64_t end = file->tell ();

            // we extrapolate the result
            number    = (number    * end) / current;
            maxSize   = (maxSize   * end) / current;
            totalSize = (totalSize * end) / current;
        }

        // we clean up resources
        delete file;
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
