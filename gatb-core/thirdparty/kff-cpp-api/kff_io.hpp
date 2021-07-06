/**
 * @file kff_io.cpp
 *
 * @brief This file define the classes needed to create and read kff files.
 * This file contains a low level and a high level APIs.
 * The Kff_file class is the base class for the low level API.
 *
 * @author Yoann Dufresne
 * Contact: yoann.dufresne0@gmail.com
 *
 */

#include <fstream>
#include <unordered_map>
#include <vector>

#ifndef KFF_IO
#define KFF_IO

#ifdef _WIN32
#include <iso646.h>
using uint = unsigned long;
#endif

// the configured options and settings for Tutorial
#define KFF_VERSION_MAJOR 1
#define KFF_VERSION_MINOR 0

class Section_GV;
class Section_Raw;
class Section_Minimizer;
class Section_Index;
class Kff_reader;

/**
 * This class is the central class for the low level kff file API.
 *
 * Kff_file directly contains the functions needed to create/read the kff file headers.
 * The class also contains functions to open section objects (one per section type).
 *
 */
class Kff_file {
private:
	bool is_writer;
	bool is_reader;
	bool tmp_closed;
	std::string filename;
	long end_position;

	// True when the header have been read/write and the fs pointer is after.
	bool header_over;

	friend class Section;
	friend class Section_GV;
	friend class Section_Index;
	friend class Section_Raw;
	friend class Section_Minimizer;
	friend class Kff_reader;

	/**
	 * Read encoding from file and save it to the public argument "encoding".
	 */
	void read_encoding();

	/**
	 * Read the metadata size from the file and save the value in metadate_size.
	 */
	void read_size_metadata();

	bool footer_discovery_ended;
	bool index_discovery_ended;

	void write_footer();
	void footer_discovery();
	void index_discovery();
	void read_index(long position);

public:
	uint8_t major_version;
	uint8_t minor_version;
	bool uniqueness;
	bool canonicity;

	Section_GV * footer;
	std::vector<Section_Index *> index;

	bool indexed;
	std::map<int64_t, char> section_positions;

	std::fstream fs;
	// encoding:        A:0  C:1 G:3 T:2
	uint8_t encoding[4] = {0, 1, 3, 2};
	
	uint32_t metadata_size = 0;

	std::unordered_map<std::string, uint64_t> global_vars;

	// --- General functions ---
	/** Open the file filename with the mode mode.
	 * mode must be choose in the set of values {r: read, w: write}
	 *
	 * @param filename The path to the file to construct/read.
	 * @param mode Opening mode of the file. w for writing, r for reading.
   *
	 */
	Kff_file(const std::string filename, const std::string mode);
	/**
	 * Destroy properly a Kff object, closing the file behind it.
	 */
	~Kff_file();
	/**
	 * In writing mode, the KFF files are indexed by default.
	 * This function allows the deactivation of the footer index generation.
	 */
	void set_indexation(bool indexed);
	/**
	 * Register a section into index
	 */
	void register_position(char section_type);
	/**
	 * Release the file pointer by temporarily close the file stream.
	 * The usage of this function increase the execution time.
	 * It should be use only when the user need a large amount of simulataneous kff files.
	 */
	void tmp_close();
	/**
	 * Reopen the file and point to the end of it.
	 * This function is automatically called on write events if the file has been temporarily closed.
	 */
	void reopen();
	/** 
	 * Close the file
	 *
	 * Physically close the file and save buffered parts to the disk if w.
	 */
	void close();

	// --- header functions ---
	/**
	 * Set the encoding used to compact the nucleotides into 2-bits values.
	 * Only the two lower bits of each uint8_t will be used.
	 * The 4 2-bits values must be diferent to each other.
	 *
	 * @param a encoding for A or a letters.
	 * @param c encoding for C or c letters.
	 * @param g encoding for G or g letters.
	 * @param t encoding for T or t letters.
	 */
	void write_encoding(uint8_t a, uint8_t c, uint8_t g, uint8_t t);
	/**
	 * Set the encoding used to compact the nucleotides into 2-bits values.
	 * Only the two lower bits of each uint8_t will be used.
	 * The 4 2-bits values must be diferent to each other.
	**/
	void write_encoding(uint8_t * encoding);
	/** Set the uniqueness of the kmers in the file.
	  * true means that no kmer will be present 2 times in the file.
	  **/
	void set_uniqueness(bool uniqueness);
	/** Set the cononicity of the kmers in the file.
	  * true means that if a kmer is present in the file, its reverse complement is absent.
	  **/
	void set_canonicity(bool canonicity);
	/**
	 * Write the metadata entered by the user.
	 * Also write all the needed values like metadata section size into the file.
	 *
	 * @param size The size of data array.
	 * @param data An array of data. Can either be plain text or binary data.
	 *
	 */
	void write_metadata(uint32_t size, const uint8_t * data);
	/**
	 * Read the next size Bytes in the file and return them into the data array.
	 * This function is used to read the metadata field.
	 * The size of the metadata is saved in the metadata_size field.
	 *
	 * @param data Array filled with the file content (up to size).
	 *
	 */
	void read_metadata(uint8_t * data);
	/**
	 * Function called to complete a header reading or writing before the first section io
	 */
	void complete_header();


	// --- general section ---
	/**
	 * Read the next Byte and return it.
	 * If the file pointer is correctly aligned with the beginning of a section, this Byte is a char that represents the section type.
	 *
	 * @return A char corresponding to the type of the following section.
	 *
	 */
	char read_section_type();
	/**
	 * Jump over the next section if it is section containing kmers.
	 *
	 * @return True if a section is skipped, false if the section can't be skipped, eof is reached or the file is not read ready.
	 */
	bool jump_next_section();
};


class Section {
	public:
		Kff_file * file;
		long beginning;

		Section(Kff_file * file);
		virtual ~Section() {};
		void close();
};

/**
 * File manipulator for Global Variable sections.
 *
 */
class Section_GV : public Section {
private:
	void read_section();
	void read_var();

	// friend class Kff_file;

public:
	uint64_t nb_vars;
	/**
	 * A map containing all the declared variables for this section.
	 * If the file is opened in r mode, all the variables of this section have been loaded during the object construction.
	 */
	std::map<std::string, uint64_t> vars;

	Section_GV(Kff_file * file);
	/**
	 * Write a variable into the file.
	 *
	 * @param var_name String name of the variable.
	 * @param value The 64 bits value of the constant.
	 *
	 */
	void write_var(const std::string & var_name, uint64_t value);
	/**
	 * Close the section.
	 * If w mode, go back to the beginning of the section to write the correct number of variables.
	 *
	 */
	void close();
};

class Block_section_reader {
protected:
	uint64_t k;
	uint64_t max;
	uint64_t data_size;
	uint8_t nb_kmers_bytes;

	friend class Kff_reader;

public:
	/**
	 * The number of blocks in the section.
	 * mode r: filled during the construction of the object.
	 * mode w: Increased each time a block is added. Used on close to write the number at the beginning of the section.
	 */
	uint64_t nb_blocks;
	uint64_t remaining_blocks;

	virtual ~Block_section_reader(){};

	static Block_section_reader * construct_section(Kff_file * file);
	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq array and the data in the data array.
	 * These arrays must be already allocated.
	 * 
	 * @param seq Array filled with the compacted sequence (2 bit / nucl).
	 * @param data Array filled with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 *
	 */
	virtual uint64_t read_compacted_sequence(uint8_t* seq, uint8_t* data) {return 0;}
	/**
	 * Jumb over the next block of the section.
	 */
	virtual void jump_sequence() {}
	/**
	 * Jump over the full section
	 */
	void jump_section() {
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}
};

/**
 * File manipulator for Raw sections.
 *
 */
class Section_Index : public Section {
private:

	friend class Kff_file;

public:
	std::map<int64_t, char> index;
	int64_t next_index;

	Section_Index(Kff_file * file);
	~Section_Index(){};
	void register_section(char section_type, int64_t index);
	void set_next_index(int64_t index);
	void close();
};

/**
 * File manipulator for Raw sections.
 *
 */
class Section_Raw: public Section, public Block_section_reader {
private:
	friend class Kff_file;
	friend class Block_section_reader;

	uint32_t read_section_header();

public:
	Section_Raw(Kff_file * file);
	~Section_Raw(){};

	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq array and the data in the data array.
	 * These arrays must be already allocated.
	 * 
	 * @param seq Array filled with the compacted sequence (2 bit / nucl).
	 * @param data Array filled with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 *
	 */
	uint64_t read_compacted_sequence(uint8_t* seq, uint8_t* data);
	/**
	 * Write a block containing the sequence seq with kmer data associated.
	 * 
	 * @param seq A compacted sequence (2 bit / nucl).
	 * @param seq_size Size of the sequence (in nucleotides).
	 * @param data Data array of the kmers in the sequence.
	 *
	 */
	void write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint8_t * data_array);
	/**
	 * Jumb over the next block of the section.
	 */
	void jump_sequence();
	/**
	 * Close the section.
	 * If w mode, go back to the beginning of the section to write the correct number of blocks.
	 *
	 */
	void close();
};


/**
 * File manipulator for Minimizer sections.
 *
 */
class Section_Minimizer : public Section, public Block_section_reader {
private:
	Section_Minimizer(Section_Minimizer & sm) = delete;
	Section_Minimizer() = delete;

	friend class Kff_file;
	friend class Block_section_reader;

	uint8_t m;
	uint8_t nb_bytes_mini;
	uint8_t mini_pos_bytes;

	uint32_t read_section_header();

public:
	Section_Minimizer(Kff_file * file);
	Section_Minimizer& operator= ( Section_Minimizer && sm);
	~Section_Minimizer();
	/**
	 * Compacted minimizer
	 */
	uint8_t * minimizer;

	/**
	 * Write the minimizer at the beginning of the block.
	 * 
	 * @param minimizer Compacted minimizer (2 bits/nucleotide)
	 *
	 */
	void write_minimizer(uint8_t * minimizer);
	/**
	 * Write a block where the sequence is composed of the nucleotides that are not part of the minimizer.
	 *
	 * @param seq Sequence without the minimizer nucleotides.
	 * @param seq_size The sequence size without its minimizer.
	 * @param mini_pos Position where the minimizer should be inserted in the sequence.
	 * @param data_array Array filled with data linked to the kmers of the sequence.
	 *
	 */
	void write_compacted_sequence_without_mini(uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array);
	/**
	 * Write a block for the couple seq/data_array.
	 * Before writing, the sequence is processed to remove the minimizer from it.
	 *
	 * @param seq Sequence without the minimizer nucleotides.
	 * @param seq_size The sequence size without its minimizer.
	 * @param mini_pos Start position of the minimizer.
	 * @param data_array Array filled with data linked to the kmers of the sequence.
	 *
	 */
	void write_compacted_sequence (uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array);
	/**
		* Add the minimizer of the block indide of a sequence.
		*
		* @param nb_kmer Number of kmers inside of the sequence
		* @param seq Compacted sequence without minimizer. The space allocaded must be
		* sufficient for the minimizer add.
		* @param mini_pos the position of the minimizer inside of the sequence.
	 	**/
	void add_minimizer(uint64_t nb_kmer, uint8_t * seq, uint64_t mini_pos);
	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq array and the data in the data array.
	 * These arrays must be already allocated.
	 * 
	 * @param seq Array filled with the compacted sequence (2 bit / nucl) where the minimizer is absent.
	 * @param data Array filled with data linked to the kmers of the sequence.
	 * @param mini_pos Index of the first character of the minimizer.
	 *
	 * @return The number of kmers in the sequence.
	 *
	 */
	uint64_t read_compacted_sequence_without_mini(uint8_t* seq, uint8_t* data, uint64_t & mini_pos);
	/**
	 * Read the next block of the section.
	 * The sequence of the block is pushed in the seq array and the data in the data array.
	 * These arrays must be already allocated.
	 * 
	 * @param seq Array filled with the compacted sequence (2 bit / nucl).
	 * @param data Array filled with data linked to the kmers of the sequence.
	 *
	 * @return The number of kmers in the sequence.
	 *
	 */
	uint64_t read_compacted_sequence(uint8_t* seq, uint8_t* data);
	/**
	 * Jumb over the next block of the section.
	 */
	void jump_sequence();
	/**
	 * Close the section.
	 * If w mode, go back to the beginning of the section to write the correct number of blocks.
	 *
	 */
	void close();
};


class Kff_reader {
private:
	Kff_file * file;
	// Space alocated for copying the current kmer given to the user.
	uint8_t * current_kmer;
	// Current sequence
	uint8_t * current_sequence;
	// Current sequence shifted to match the 4 different alignements
	uint8_t ** current_shifts;
	// Size in nucleotides of current sequence
	uint64_t current_seq_nucleotides;
	// Size in bytes of current sequence
	uint64_t current_seq_bytes;
	// Number of kmers in the current sequence
	uint64_t current_seq_kmers;
	// Number of kmer remaining in the current block
	uint64_t remaining_kmers;
	// Data array for the current block
	uint8_t * current_data;
	// Section currently readed.
	Block_section_reader * current_section;
	// Remaining blocks before end of the section
	uint64_t remaining_blocks;

	void read_until_first_section_block();
	void read_next_block();

public:
	
	Kff_reader(std::string filename);
	~Kff_reader();

	bool has_next();
	uint64_t next_block(uint8_t* & sequence, uint8_t* & data);
	void next_kmer(uint8_t* & sequence, uint8_t* & data);

	uint64_t get_var(std::string name);
	uint8_t * get_encoding();
};

#endif
