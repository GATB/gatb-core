#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>
#include <math.h>

#include <map>

#include "kff_io.hpp"

using namespace std;


// Utils

template<typename T>
void write_value(T val, fstream & fs) {
	fs.write((char*)&val, sizeof(T));
}
template<typename T>
void read_value(T & val, fstream & fs) {
	fs.read((char*)&val, sizeof(T));
}

uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem) {
	if (bits_per_elem == 0 or nb_elem == 0)
		return 0;
	else
		return ((bits_per_elem * nb_elem - 1) / 8) + 1;
}


// ----- Open / Close functions -----

Kff_file::Kff_file(const string filename, const string mode) {
	// Variable init
	this->is_writer = false;
	this->is_reader = false;
	auto streammode = fstream::binary;

	// Determine the mode and open the file
	if (mode[0] == 'w') {
		this->is_writer = true;
		streammode |= fstream::out;
	} else if (mode[0] == 'r') {
		this->is_reader = true;
		streammode |= fstream::in;		
	} else {
		cerr << "Unsupported mode " << mode << endl;
		exit(0);
	}

	// Open the file
	this->filename = filename;
	this->fs.open(filename, streammode);
	this->tmp_closed = false;
	this->header_over = false;

	// Write the version at the beginning of the file
	if (this->is_writer) {
		this->fs << (char)KFF_VERSION_MAJOR << (char)KFF_VERSION_MINOR;
	}

	// Read the header
	else if (this->is_reader) {
		this->fs >> this->major_version >> this->minor_version;

		if (KFF_VERSION_MAJOR < this->major_version or (KFF_VERSION_MAJOR == this->major_version and KFF_VERSION_MINOR < this->minor_version)) {
			cerr << "The software version " << (uint)KFF_VERSION_MAJOR << "." << (uint)KFF_VERSION_MINOR << " can't read files writen in version " << (uint)this->major_version << "." << (uint)this->minor_version << endl;
			throw "Unexpected version number";
		}

		this->read_encoding();
		this->read_size_metadata();
	}
}


Kff_file::~Kff_file() {
	this->close();
}


void Kff_file::complete_header() {
	if (this->header_over)
		return;

	// If the metadata has not been read, jump over
	if (this->is_reader) {
		this->fs.seekp(this->fs.tellp() + (long)this->metadata_size);
	}

	// If metadata has not been write, write a 0 byte one.
	else if (this->is_writer) {
		this->write_metadata(0, nullptr);
	}

	this->header_over = true;
}


void Kff_file::tmp_close() {
	if (this->is_writer and this->fs.is_open()) {
		this->fs.close();
		this->fs.clear();
		this->tmp_closed = true;
	}
}


void Kff_file::reopen() {
	if (this->tmp_closed) {
		auto streammode = fstream::binary | fstream::out | fstream::in | fstream::ate;

		// Open the file
		this->fs.open(this->filename, streammode);
		this->tmp_closed = false;
	}
}


void Kff_file::close() {
	if (this->fs.is_open())
		this->fs.close();
	this->is_writer = false;
	this->is_reader = false;
}


// ----- Header functions -----

void Kff_file::write_encoding(uint8_t a, uint8_t c, uint8_t g, uint8_t t) {
	// Value masking
	a &= 0b11; c &= 0b11; g &= 0b11; t &= 0b11;

	// Verification of the differences.
	assert(a != c); assert(a != g); assert(a != t);
	assert(c != g); assert(g != t);
	assert(g != t);

	// set values
	this->encoding[0] = a;
	this->encoding[1] = c;
	this->encoding[2] = g;
	this->encoding[3] = t;

	// Write to file
	uint8_t code = (a << 6) | (c << 4) | (g << 2) | t;
	this->fs << code;
}

void Kff_file::write_encoding(uint8_t * encoding) {
	this->write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
}

void Kff_file::read_encoding() {
	uint8_t code, a, c, g, t;
	// Get code values
	this->fs >> code;

	// Split each nucleotide encoding
	this->encoding[0] = a = (code >> 6) & 0b11;
	this->encoding[1] = c = (code >> 4) & 0b11;
	this->encoding[2] = g = (code >> 2) & 0b11;
	this->encoding[3] = t = code & 0b11;

	// Verification of the encoding
	if (a == c or a == g or a == t or c == g or c == t or g == t) {
		throw "Wrong encoding. The 4 2-bits values must be different.";
	}
}

void Kff_file::write_metadata(uint32_t size, const uint8_t * data) {
	write_value(size, fs);
	this->fs.write((char *)data, size);
	this->header_over = true;
}

void Kff_file::read_size_metadata() {
	read_value(this->metadata_size, fs);
}

void Kff_file::read_metadata(uint8_t * data) {
	this->fs.read((char *)data, this->metadata_size);
	this->header_over = true;
}

bool Kff_file::jump_next_section() {
	if (not is_reader)
		return false;
	char section_type = read_section_type();
	if (fs.eof())
		return false;
	if (section_type == 'r' or section_type == 'm') {
		Block_section_reader * section = Block_section_reader::construct_section(this);
		section->jump_section();
		delete section;
		return true;
	}
	return false;
}


// ----- Sections -----

char Kff_file::read_section_type() {
	// Verify that header has been read.
	if (not this->header_over) {
		this->complete_header();
	}

	char type = '\0';
	this->fs >> type;
	this->fs.seekp(this->fs.tellp() - 1l);
	return type;
}


// ----- Global variables sections -----

Section_GV::Section_GV(Kff_file * file) {
	if (file->tmp_closed) {
		file->reopen();
	}

	if (not file->header_over) {
		file->complete_header();
	}

	this->file = file;
	this->begining = file->fs.tellp();
	this->nb_vars = 0;
	this->is_closed = true;

	if (file->is_reader) {
		this->read_section();
	}

	if (file->is_writer) {
		this->is_closed = false;
		fstream & fs = file->fs;
		fs << 'v';
		write_value(nb_vars, fs);
	}
}

void Section_GV::write_var(const string & var_name, uint64_t value) {
	assert(!this->is_closed);

	if (file->tmp_closed) {
		file->reopen();
	}

	auto & fs = this->file->fs;
	fs << var_name << '\0';
	write_value(value, fs);
	this->nb_vars += 1;

	this->vars[var_name] = value;
	this->file->global_vars[var_name] = value;
}

void Section_GV::read_section() {
	char type;
	file->fs >> type;
	if (type != 'v')
		throw "The section do not start with the 'v' char, you can't open a Global Variable section.";

	read_value(this->nb_vars, file->fs);
	for (uint64_t i=0 ; i<nb_vars ; i++) {
		this->read_var();
	}
}

void Section_GV::read_var() {
	if (file->fs.eof())
		throw "eof reached before the end of the variable section";

	// Name reading
	stringstream ss;
	char c = 'o';
	this->file->fs >> c;
	while (c != '\0') {
		ss << c;
		this->file->fs >> c;
	}

	// Value reading
	uint64_t value;
	read_value(value, file->fs);

	// Saving
	string name = ss.str();
	this->vars[name] = value;
	this->file->global_vars[name] = value;
}

void Section_GV::close() {
	if (file->is_writer) {
		if (file->tmp_closed)
			file->reopen();
		// Save current position
		fstream &	 fs = this->file->fs;
		long position = fs.tellp();
		// Go write the number of variables in the correct place
		fs.seekp(this->begining + 1);
		write_value(nb_vars, fs);

		fs.seekp(position);
		this->is_closed = true;
	}
}

Block_section_reader * Block_section_reader::construct_section(Kff_file * file) {
	// Very and complete if needed the header
	file->complete_header();

	char type = file->read_section_type();
	if (type == 'r') {
		return new Section_Raw(file);
	} else if (type == 'm') {
		return new Section_Minimizer(file);
	} else
		return nullptr;
}


// ----- Raw sequence section -----

Section_Raw::Section_Raw(Kff_file * file) {
	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing k variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing data_size variable";
	
	uint64_t k = file->global_vars["k"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	if (file->tmp_closed) {
		file->reopen();
	}

	if (not file->header_over) {
		file->complete_header();
	}

	this->file = file;
	this->begining = file->fs.tellp();
	this->nb_blocks = 0;
	this->is_closed = true;

	this->k = k;
	this->max = max;
	this->data_size = data_size;

	// Computes the number of bytes needed to store the number of kmers in each block
	uint64_t nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = static_cast<uint8_t>(bytes_from_bit_array(nb_bits, 1));

	if (file->is_reader) {
		this->read_section_header();
	}

	if (file->is_writer) {
		this->is_closed = false;
		fstream & fs = file->fs;
		fs << 'r';
		write_value(nb_blocks, fs);
	}
}

uint32_t Section_Raw::read_section_header() {
	fstream & fs = file->fs;

	char type;
	fs >> type;
	if (type != 'r')
		throw "The section do not start with the 'r' char, you can't open a Raw sequence section.";

	read_value(this->nb_blocks, fs);
	this->remaining_blocks = this->nb_blocks;

	return this->nb_blocks;
}

void Section_Raw::write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint8_t * data_array) {
	assert(!this->is_closed);
	if (file->tmp_closed) {
		file->reopen();
	}
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size - k + 1;
	file->fs.write((char*)&nb_kmers, this->nb_kmers_bytes);
	// 2 - Write sequence
	uint64_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	this->file->fs.write((char *)seq, seq_bytes_needed);
	// 3 - Write data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers);
	this->file->fs.write((char *)data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Raw::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0)
		file->fs.read((char*)&nb_kmers_in_block, this->nb_kmers_bytes);
	// 2 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	file->fs.read((char*)seq, seq_bytes_needed);
	// 3 - Read the data.
	uint64_t data_bytes_used = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	file->fs.read((char*)data, data_bytes_used);

	this->remaining_blocks -= 1;

	return nb_kmers_in_block;
}


void Section_Raw::jump_sequence() {
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0)
		file->fs.read((char*)&nb_kmers_in_block, this->nb_kmers_bytes);
	// 2 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	// 3 - Determine the data size
	size_t data_bytes_used = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	// 4 - Jumb over the 
	file->fs.seekp(file->fs.tellp() + (long)(seq_bytes_needed + data_bytes_used));
	this->remaining_blocks -= 1;
}


void Section_Raw::close() {
	if (file->is_writer) {
		if (file->tmp_closed) {
			file->reopen();
		}
		// Save current position
		fstream &	 fs = this->file->fs;
		long position = fs.tellp();
		// Go write the number of variables in the correct place
		fs.seekp(this->begining + 1);
		write_value(nb_blocks, fs);
		fs.seekp(position);
		this->is_closed = true;
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}
}



// ----- Minimizer sequence section -----

Section_Minimizer::Section_Minimizer(Kff_file * file) {
	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing k variable";
	if (file->global_vars.find("m") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing m variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the minimizer section due to missing data_size variable";
	
	uint64_t k = file->global_vars["k"];
	uint64_t m = file->global_vars["m"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	if (file->tmp_closed) {
		file->reopen();
	}

	if (not file->header_over) {
		file->complete_header();
	}

	this->file = file;
	this->begining = file->fs.tellp();
	this->nb_blocks = 0;
	this->is_closed = true;

	this->k = k;
	this->m = m;
	this->max = max;
	this->data_size = data_size;

	// Computes the number of bytes needed to store the number of kmers in each block
	uint64_t nb_bits = static_cast<uint64_t>(ceil(log2(max)));
	this->nb_kmers_bytes = static_cast<uint8_t>(bytes_from_bit_array(nb_bits, 1));
	this->nb_bytes_mini = static_cast<uint8_t>(bytes_from_bit_array(2, m));
	this->minimizer = new uint8_t[nb_bytes_mini];
	memset(this->minimizer, 0, nb_bytes_mini);
	uint64_t mini_pos_bits = static_cast<uint8_t>(ceil(log2(k+max-1)));
	this->mini_pos_bytes = bytes_from_bit_array(mini_pos_bits, 1);

	if (file->is_reader) {
		this->read_section_header();
	}

	if (file->is_writer) {
		this->is_closed = false;
		fstream & fs = file->fs;
		fs << 'm';
		this->write_minimizer(this->minimizer);
		file->fs.seekp(file->fs.tellp()+(long)nb_bytes_mini);
		write_value(nb_blocks, fs);
	}
}

Section_Minimizer::~Section_Minimizer() {
	delete[] this->minimizer;
}

Section_Minimizer& Section_Minimizer::operator= ( Section_Minimizer && sm) {
	file = sm.file;
	sm.file = nullptr;
	begining = sm.begining;
	is_closed = sm.is_closed;
	nb_blocks = sm.nb_blocks;
	m = sm.m;
	nb_bytes_mini = sm.nb_bytes_mini;
	std::swap(minimizer, sm.minimizer);

	return *this;
}

void Section_Minimizer::write_minimizer(uint8_t * minimizer) {
	assert(!this->is_closed);
	if (file->tmp_closed) {
		file->reopen();
	}

	uint64_t pos = file->fs.tellp();
	file->fs.seekp(this->begining+1);
	file->fs.write((char *)minimizer, this->nb_bytes_mini);
	memcpy(this->minimizer, minimizer, this->nb_bytes_mini);
	file->fs.seekp(pos);
	
}

uint32_t Section_Minimizer::read_section_header() {
	fstream & fs = file->fs;

	// Verify section type
	char type;
	fs >> type;
	if (type != 'm')
		throw "The section do not start with the 'm' char, you can't open a Minimizer sequence section.";

	// Read the minimizer
	file->fs.read((char *)this->minimizer, this->nb_bytes_mini);

	// Read the number of following blocks
	read_value(nb_blocks, fs);
	this->remaining_blocks = this->nb_blocks;
	return nb_blocks;
}

void Section_Minimizer::write_compacted_sequence_without_mini(uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array) {
	assert(!this->is_closed);
	if (file->tmp_closed) {
		file->reopen();
	}
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size + m - k + 1;
	file->fs.write((char*)&nb_kmers, this->nb_kmers_bytes);
	// 2 - Write minimizer position
	file->fs.write((char *)&mini_pos, this->mini_pos_bytes);
	// 3 - Write sequence with chopped minimizer
	uint64_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	this->file->fs.write((char *)seq, seq_bytes_needed);
	// 4 - Write data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers);
	this->file->fs.write((char *)data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Minimizer::read_compacted_sequence_without_mini(uint8_t* seq, uint8_t* data, uint64_t & mini_pos) {
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0)
		file->fs.read((char*)&nb_kmers_in_block, this->nb_kmers_bytes);
	// 2 - Read the minimizer position
	uint64_t tmp_mini_pos = 0;
	file->fs.read((char *)&tmp_mini_pos, this->mini_pos_bytes);
	mini_pos = tmp_mini_pos;
	// 3 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	file->fs.read((char*)seq, seq_bytes_needed);
	// 4 - Read the data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	file->fs.read((char*)data, data_bytes_needed);
	// cout << data_bytes_needed << endl;

	this->remaining_blocks -= 1;
	return nb_kmers_in_block;
}

void Section_Minimizer::jump_sequence() {
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0)
		file->fs.read((char*)&nb_kmers_in_block, this->nb_kmers_bytes);
	// 2 - Read the minimizer position
	uint64_t tmp_mini_pos = 0;
	file->fs.read((char *)&tmp_mini_pos, this->mini_pos_bytes);
	// 3 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	// 3 - Determine the data size
	size_t data_bytes_used = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	// 4 - Jumb over the 
	file->fs.seekp(file->fs.tellp() + (long)(seq_bytes_needed + data_bytes_used));

	this->remaining_blocks -= 1;
}


/* Bitshift to the left all the bits in the array with a maximum of 7 bits.
 * Overflow on the left will be set into the previous cell.
 */
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=0 ; i<length-1 ; i++) {
		bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
	}
	bitarray[length-1] <<= bitshift;
}

/* Similar to the previous function but on the right */
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=length-1 ; i>0 ; i--) {
		bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
	}
	bitarray[0] >>= bitshift;
}

/* Fusion to bytes into one.
 * The merge_index higher bits are from left_bits the others from right_bits
 */
inline uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}

void Section_Minimizer::write_compacted_sequence (uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array) {
	assert(!this->is_closed);
	if (file->tmp_closed) {
		file->reopen();
	}
	// Compute all the bit and byte quantities needed.
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	uint left_offset_nucl = (4 - seq_size % 4) % 4;

	// 1 - Prepare space for sequence manipulation
	uint8_t * seq_copy = new uint8_t[seq_bytes];
	memcpy(seq_copy, seq, seq_bytes);
	
	// 2 - Move the suffix to the new bytes
	uint mini_start_byte = (mini_pos + left_offset_nucl) / 4;
	uint suff_start_byte = (mini_pos + m + left_offset_nucl) / 4;
	for (uint i=0 ; i<seq_bytes-suff_start_byte ; i++) {
		seq_copy[mini_start_byte + i] = seq_copy[suff_start_byte + i];
	}

	// 3 - shift the suffix to the right bit
	uint mini_offset = (mini_pos + left_offset_nucl) % 4;
	uint suff_offset = (mini_pos + m + left_offset_nucl) % 4;
	if (mini_offset < suff_offset)
		leftshift8(seq_copy + mini_start_byte, seq_bytes - mini_start_byte, (suff_offset - mini_offset) * 2);
	else
		rightshift8(seq_copy + mini_start_byte, seq_bytes - mini_start_byte, (mini_offset - suff_offset) * 2);

	// 4 - fusion the common byte
	seq_copy[mini_start_byte] = fusion8(seq[mini_start_byte], seq_copy[mini_start_byte], mini_offset * 2);

	// 5 - put all the offset bits on the left
	leftshift8(seq_copy, seq_bytes, left_offset_nucl * 2);
	rightshift8(seq_copy, seq_bytes, ((4 - ((seq_size - m) % 4)) % 4) * 2);

	// Write the compacted sequence into file.
	this->write_compacted_sequence_without_mini(seq_copy, seq_size-m, mini_pos, data_array);

	delete[] seq_copy;
}

void Section_Minimizer::add_minimizer(uint64_t nb_kmer, uint8_t * seq, uint64_t mini_pos) {
	uint64_t seq_size = nb_kmer + k - 1;
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	uint64_t size_no_mini = seq_size - m;
	uint64_t bytes_no_mini = bytes_from_bit_array(2, size_no_mini);
	uint64_t left_offset = (4 - size_no_mini) % 4;
	// Shift the whole sequence to the left to have np padding on byte 0.
	leftshift8(seq, bytes_no_mini, left_offset*2);

	// Compute usefull quantities
	uint64_t nucl_mini_start = mini_pos;
	uint64_t byte_mini_start = nucl_mini_start / 4;
	uint8_t offset_mini_start = 2 * (nucl_mini_start % 4);

	uint64_t nucl_mini_stop = nucl_mini_start + m - 1;
	uint8_t offset_mini_stop = 2 * (nucl_mini_stop % 4);

	// Copy the suffix to shift it
	uint64_t size_suffix = size_no_mini - mini_pos;
	uint64_t suffix_bytes = bytes_from_bit_array(2, size_suffix);
	uint64_t suffix_start_byte = byte_mini_start;
	uint64_t suffix_stop_nucl = nucl_mini_start + size_suffix - 1;
	uint64_t suffix_stop_byte = suffix_stop_nucl / 4;

	uint8_t * suffix = new uint8_t[suffix_bytes + 1];
	memset(suffix, 0, suffix_bytes + 1);
	memcpy(suffix, seq+suffix_start_byte, suffix_stop_byte - suffix_start_byte + 1);
	leftshift8(suffix, suffix_stop_byte - suffix_start_byte + 1, offset_mini_start);
	rightshift8(suffix, suffix_stop_byte - suffix_start_byte + 1, (offset_mini_stop+2)%8);

	// // Add the suffix
	uint64_t suff_first_nucl_used = offset_mini_stop/2 + 1;
	uint64_t suff_first_byte_used = suff_first_nucl_used / 4;
	uint64_t suff_last_nucl_used = suff_first_nucl_used + size_suffix - 1;
	uint64_t suff_last_byte_used = suff_last_nucl_used / 4;
	uint64_t bytes_suff_used = suff_last_byte_used - suff_first_byte_used + 1;
	
	uint64_t fusion_byte = (nucl_mini_stop + 1) / 4;
	seq[fusion_byte] = fusion8(seq[fusion_byte], suffix[0], 2 * (suff_first_nucl_used % 4));
	for (uint64_t i=1 ; i<=bytes_suff_used ; i++) {
		seq[fusion_byte+i] = suffix[i];
	}

	// Prepare the minimizer
	uint8_t * mini_shifted = new uint8_t[this->nb_bytes_mini+1];
	mini_shifted[this->nb_bytes_mini] = 0;
	memcpy(mini_shifted, this->minimizer, this->nb_bytes_mini);
	uint64_t mini_left_offset = (8 - ((2 * m) % 8)) % 8;

	leftshift8(mini_shifted, this->nb_bytes_mini, mini_left_offset);
	rightshift8(mini_shifted, this->nb_bytes_mini+1, offset_mini_start);
	// Add the minimizer inside the sequence
	seq[byte_mini_start] = fusion8(seq[byte_mini_start], mini_shifted[0], offset_mini_start);
	for (uint64_t i=1 ; i<fusion_byte-byte_mini_start ; i++) {
		seq[byte_mini_start+i] = mini_shifted[i];
	}

	seq[fusion_byte] = fusion8(mini_shifted[fusion_byte-byte_mini_start], seq[fusion_byte], 2 * (suff_first_nucl_used % 4));

	rightshift8(seq, seq_bytes, (8 - 2 * (seq_size % 4)) % 8);

	delete[] suffix;
	delete[] mini_shifted;
}

uint64_t Section_Minimizer::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	// Read the block
	uint64_t mini_pos;
	uint64_t nb_kmers_in_block = this->read_compacted_sequence_without_mini(seq, data, mini_pos);
	this->add_minimizer(nb_kmers_in_block, seq, mini_pos);
	
	return nb_kmers_in_block;
}

void Section_Minimizer::close() {
	if (file->is_writer) {
		if (file->tmp_closed) {
			file->reopen();
		}
		// Save current position
		fstream &	fs = this->file->fs;
		long position = fs.tellp();
		// Go write the number of variables in the correct place
		fs.seekp(this->begining + 1l + (long)this->nb_bytes_mini);
		write_value(nb_blocks, fs);
		fs.seekp(position);
		this->is_closed = true;
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}
}





// -------- Start of the high level API -----------

Kff_reader::Kff_reader(std::string filename) {
	// Open the file
	this->file = new Kff_file(filename, "r");

	// Create fake small datastrucutes waiting for the right values.
	this->current_shifts = new uint8_t*[4];
	for (uint8_t i=0 ; i<4 ; i++) {
		this->current_shifts[i] = new uint8_t[1];
		this->current_shifts[i][0] = 0;
	}
	this->current_sequence = this->current_shifts[0];
	this->current_data = new uint8_t[1];
	this->current_data[0] = 0;

	this->current_section = NULL;
	this->current_kmer = new uint8_t[1];
	this->remaining_kmers = 0;

	this->has_next();
}

Kff_reader::~Kff_reader() {
	delete[] this->current_kmer;
	delete[] this->current_data;

	for (uint i=0 ; i<4 ; i++)
		delete[] this->current_shifts[i];
	delete[] this->current_shifts;

	delete this->file;
}

void Kff_reader::read_until_first_section_block() {
	while (current_section == NULL or remaining_blocks == 0) {
		// char section_type = this->file->read_section_type();
		char section_type = file->read_section_type();
		if (section_type == 0)
			if (this->file->fs.eof())
				break;
		// --- Update data structure sizes ---
		if (section_type == 'v') {
			// Read the global variable block
			Section_GV gvs(file);
			// Update sequence size if k or max change
			if (gvs.vars.find("k") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end()) {
				// Compute the max size of a sequence
				auto k = this->file->global_vars["k"];
				auto max = this->file->global_vars["max"];
				uint64_t max_size = bytes_from_bit_array(2, max + k - 1);
				// Allocate the right amount of memory and place the pointers to the right addresses
				for (uint8_t i=0 ; i<4 ; i++) {
					delete[] this->current_shifts[i];
					this->current_shifts[i] = new uint8_t[max_size];
					memset(this->current_shifts[i], 0, max_size);
				}
				this->current_sequence = this->current_shifts[0];
				delete[] this->current_kmer;
				this->current_kmer = new uint8_t[k/4 + 1];
				memset(this->current_kmer, 0, (k/4+1));
			}

			// Update the data array size
			if (gvs.vars.find("data_size") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end()) {
				// Compute the max size of a data array
				auto data_size = this->file->global_vars["data_size"];
				auto max = this->file->global_vars["max"];
				uint64_t max_size = data_size * max;
				delete[] this->current_data;
				this->current_data = new uint8_t[max_size];
				memset(this->current_data, 0, max_size);
			}
		}
		// Mount data from the files to the datastructures.
		else {
			current_section = Block_section_reader::construct_section(file);
			remaining_blocks = current_section->nb_blocks;
		}
	}
}

void Kff_reader::read_next_block() {
	// Read from the file
	current_seq_kmers = remaining_kmers = current_section->read_compacted_sequence(current_sequence, current_data);
	current_seq_nucleotides = remaining_kmers + current_section->k - 1;
	current_seq_bytes = bytes_from_bit_array(2, current_seq_nucleotides);

	// Create the 4 possible shifts of the sequence for easy use.
	for (uint8_t i=1 ; i<4 ; i++) {
		// Copy
		memcpy(current_shifts[i], current_sequence, current_seq_bytes);
		// Shift
		rightshift8(current_shifts[i], current_seq_bytes, 2 * i);
	}
}

bool Kff_reader::has_next() {
	if (current_section == NULL and !file->fs.eof())
		read_until_first_section_block();
	return !file->fs.eof();
}

uint64_t Kff_reader::next_block(uint8_t* & sequence, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		sequence = NULL;
		data = NULL;
		return 0;
	}

	read_next_block();
	
	sequence = current_sequence;
	data = current_data;

	auto nb_kmers = remaining_kmers;
	remaining_kmers = 0;
	remaining_blocks -= 1;
	if (remaining_blocks == 0) {
		delete current_section;
		current_section = NULL;
	}

	return nb_kmers;
}

void Kff_reader::next_kmer(uint8_t* & kmer, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		kmer = NULL;
		data = NULL;
		return;
	}

	// Load the next block
	if (remaining_kmers == 0) {
		read_next_block();
	}

	uint64_t right_shift = (remaining_kmers - 1) % 4;
	uint64_t prefix_offset = (4 - (current_seq_nucleotides % 4)) % 4;
	uint64_t kmer_idx = current_seq_kmers - remaining_kmers;

	uint64_t start_nucl = prefix_offset + right_shift + kmer_idx;
	uint64_t start_byte = start_nucl / 4;
	uint64_t end_nucl = start_nucl + file->global_vars["k"] - 1;
	uint64_t end_byte = end_nucl / 4;

	memcpy(current_kmer, current_shifts[right_shift]+start_byte, end_byte-start_byte+1);
	kmer = current_kmer;
	data = current_data + (current_seq_kmers - remaining_kmers) * this->file->global_vars["data_size"];
	
	// Read the next block if needed.
	remaining_kmers -= 1;
	if (remaining_kmers == 0) {
		remaining_blocks -= 1;
		if (remaining_blocks == 0) {
			delete current_section;
			current_section = NULL;
		}
	}
}


uint64_t Kff_reader::get_var(string name) {
	if (file->global_vars.find(name) != file->global_vars.end())
		return file->global_vars[name];

	cerr << "Variable " << name << " is absent from the file." << endl;
	exit(2);

	return 0;
}


uint8_t * Kff_reader::get_encoding() {
	return file->encoding;
}












