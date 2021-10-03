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

#ifndef _COUNT_PROCESSOR_DUMP_KFF_HPP_
#define _COUNT_PROCESSOR_DUMP_KFF_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <kff-cpp-api/kff_io.hpp>
#include <kff-cpp-api/merge.hpp>
#include <dirent.h>
#include <errno.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** The CountProcessorDumpKff 
 *
 * see CountProcessorDump for a real doc
 *
 * */
template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorDumpKff : public CountProcessorAbstract<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type Type;
    typedef typename kmer::impl::Kmer<span>::ModelCanonical ModelCanonical;

    /** Constructor */
    CountProcessorDumpKff (
            std::string                         prefix,
        tools::storage::impl::Group&            group,
        const size_t                            kmerSize,
        system::ISynchronizer*                  synchronizer = 0,
        tools::storage::impl::Partition<Count>* solidCounts  = 0,
        size_t                                  nbPartsPerPass = 0
    )
        : _prefix(prefix), _group(group), _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass), _synchronizer(0), _model(kmerSize)
    {
        setSynchronizer (synchronizer);
        // create kff temp dir
        if(!System::file().doesExist(_prefix)){
            int ok = System::file().mkdir(_prefix, 0755);
            if(ok != 0){
                //throw Exception ("Error: can't create kff directory");
            }
        }


    }

    /** Destructor */
    virtual ~CountProcessorDumpKff ()
    {
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    void begin (const Configuration& config)
    {
        /** We remember the number of partitions for one pass. */
        _nbPartsPerPass = config._nb_partitions;

        /** We compute the number of partitions. */
        //size_t nbTotalPartitions = config._nb_partitions * config._nb_passes;
    }

    /** \copydoc ICountProcessor<span>::clones */
    CountProcessorAbstract<span>* clone ()
    {
        /** Note : we share the synchronizer for all the clones. */
        return new CountProcessorDumpKff (_prefix, _group, _kmerSize, _synchronizer, 0, _nbPartsPerPass);
    }

    /** \copydoc ICountProcessor<span>::finishClones */
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        for (size_t i=0; i<clones.size(); i++)
        {
            /** We have to recover type information. */
            if (CountProcessorDumpKff* clone = dynamic_cast<CountProcessorDumpKff*> (clones[i]))
            {
                for (std::map<std::string,size_t>::iterator it = clone->_namesOccur.begin(); it != clone->_namesOccur.end(); ++it)
                {
                    this->_namesOccur[it->first] += it->second;
                }
            }
        }

        // get list of KFF files
		std::vector<std::string> list_kff_files;
		if (auto dir = opendir((_prefix + "/").c_str())) {
			while (auto f = readdir(dir)) {
				if (!f->d_name || f->d_name[0] == '.')
					continue; // Skip everything that starts with a dot
				list_kff_files.push_back(std::string(_prefix + "/" + f->d_name));
			}
			closedir(dir);
		}
        
        // this seems like a good place to do KFF merge
        kff_merge(list_kff_files,_prefix + ".merged.kff");

        // cleanup
        for (auto kff_file : list_kff_files)
        {
            remove(kff_file.c_str());
        }

        // print list of KFF files (in case directory removal didn't work)
		if (auto dir = opendir((_prefix + "/").c_str())) {
			while (auto f = readdir(dir)) {
				if (!f->d_name || (f->d_name[0] == '.' && (f->d_name[1] == '.' || f->d_name[1] == '\0')))
					continue; // Skip "." and ".."
                std::cout << "file couldn't be deleted:" << std::string(_prefix + "/" + f->d_name) << std::endl;
			}
			closedir(dir);
		}

        int rmdir_res = rmdir((_prefix + "/").c_str());
        if (rmdir_res != 0)
            std::cout << "failed to remove folder (" << _prefix << "): "<< strerror(errno) << std::endl;
    }


    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::beginPart */
    void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name)
    {
        /** We get the actual partition idx in function of the current partition AND pass identifiers. */
        //size_t actualPartId = partId + (passId * _nbPartsPerPass);

        /** We update some stats (want to know how many "hash" or "vector" partitions we use). */
        _namesOccur[name] ++;

        _outfile = new Kff_file(_prefix + "/" +std::to_string(passId) + "-" + std::to_string(partId) + ".kff", "w");
       
        uint8_t encoding[] = {0, 1, 3, 2};
        _outfile->write_encoding(encoding);
           
        /** We save (as metadata) some information. */
        Section_GV sgv(_outfile);
        sgv.write_var("k", _kmerSize);
        sgv.write_var("max", 1); // 1 kmer per sequence
        sgv.write_var("data_size", 4); // DSK counts are stored as uint32_t
        sgv.close();

        //std::cout <<"new part " << passId << " " << partId << " " << name <<std::endl;
        sr = new Section_Raw(_outfile);
    }

    /** \copydoc ICountProcessor<span>::endPart */
    void endPart (size_t passId, size_t partId)
    {
        /** We flush the current collection for the partition just processed. */
        sr->close();
        delete _outfile;
    }

   
	// from kff-cpp-api 
		/* Transform a char * sequence into a uint8_t 2-bits/nucl
		* Encoding ACTG
		* Size must be <= 4
		*/
		uint8_t uint8_packing(std::string sequence) {
			size_t size = sequence.length();
			assert(size <= 4);

			uint8_t val = 0;
			for (size_t i=0 ; i<size ; i++) {
				val <<= 2;
				val += (sequence[i] >> 1) & 0b11;
			}

			return val;
		}

		/* Encode the sequence into an array of uint8_t packed sequence slices.
		* The encoded sequences are organised in big endian order.
		*/
		void encode_sequence(std::string sequence, uint8_t * encoded) {
			size_t size = sequence.length();
			// Encode the truncated first 8 bits sequence
			size_t remnant = size % 4;
			if (remnant > 0) {
				encoded[0] = uint8_packing(sequence.substr(0, remnant));
				encoded += 1;
			}

			// Encode all the 8 bits packed
			size_t nb_uint_needed = size / 4;
			for (size_t i=0 ; i<nb_uint_needed ; i++) {
				encoded[i] = uint8_packing(sequence.substr(remnant + 4*i, 4));
				// encoded[i] = uint8_packing(sequence + remnant + (i<<2), 4);
			}
		}


        // some questions about endianness here
        void u8from32 (uint8_t b[4], uint32_t u32)
        {
            b[3] = (uint8_t)u32;
            b[2] = (uint8_t)(u32>>=8);
            b[1] = (uint8_t)(u32>>=8);
            b[0] = (uint8_t)(u32>>=8);
        }

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
        std::string sequence = _model.toString(kmer);
		uint8_t encoded[1024];
     	uint8_t counts[4];
        encode_sequence(sequence.c_str(), encoded);
		u8from32(counts, sum);
        //std::cout << std::to_string((int)counts[0]) << "-" << std::to_string((int)counts[1]) << "-"<< std::to_string((int)counts[2])  <<"-"<< std::to_string((int)counts[3]) <<";"<<std::endl;
        sr->write_compacted_sequence(encoded, _kmerSize, counts);
        return true;
    }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getProperties */
    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;

        return result;
    }

    /** Get the number of items.
     * \return the total number of items in the partition. */
    u_int64_t getNbItems ()  { return 0; }

private:

    std::string _prefix;
    Kff_file *_outfile;
    Section_Raw *sr;


    tools::storage::impl::Group& _group;

    size_t _kmerSize;

    size_t _nbPartsPerPass;

    system::ISynchronizer* _synchronizer;
    void setSynchronizer (system::ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }

    ModelCanonical _model;

    std::map<std::string,size_t> _namesOccur;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_DUMP_HPP_ */
