#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>
#include <math.h>

#include <map>
#include <vector>

#include "kff_io.hpp"

using namespace std;


// Utils

template<typename T>
void store_big_endian(uint8_t * buff, size_t size, const T& data) {
	for (int b = size - 1; b >= 0; --b) {
		*buff++ = data >> (8 * b);
	}
}

template<typename T>
void load_big_endian(uint8_t * buff, size_t size, T& data) {
	data = 0;
	for (uint b=0 ; b < size; b++) {
		data <<= 8;
		data |= buff[b];
	}
}

uint64_t bytes_from_bit_array(uint64_t bits_per_elem, uint64_t nb_elem) {
	if (bits_per_elem == 0 or nb_elem == 0)
		return 0;
	else
		return ((bits_per_elem * nb_elem - 1) / 8) + 1;
}

static void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift);
static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift);
static uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index);


// ----- Open / Close functions -----

Kff_file::Kff_file(const string filename, const string mode) {
	// Variable init
	this->filename = filename;
	
	this->is_writer = false;
	this->is_reader = false;
	
	this->writing_started = false;
	this->next_free = 0;
	this->buffer_size = 1 << 10; // 1 KB
	// this->buffer_size = 1 << 4;
	this->max_buffer_size = 1 << 20; // 1 MB
	// this->max_buffer_size = 1 << 6;
	this->file_buffer = new uint8_t[this->buffer_size];
	this->file_size = 0;
	this->delete_on_destruction = false;

	this->open(mode);
}

void Kff_file::open(string mode) {
	this->writing_started = false;
	this->current_position = 0;

	// Determine the mode and open the file
	if (mode[0] == 'w') {
		this->is_writer = true;
		this->file_size = 0;
		this->next_free = 0;
	} else if (mode[0] == 'r') {
		this->is_reader = true;
		// If no info on the file
		if (this->file_size == 0 and this->next_free == 0) {
			// Open the fp
			this->fs.open(this->filename, fstream::binary | fstream::in);
			// Compute the file length
			long position = this->fs.tellp();
			this->fs.seekg(0, this->fs.end);
			this->file_size = (long)(this->fs.tellp()) - position;
			// Go back to the beginning
			this->fs.seekg(0, this->fs.beg);
		}
	} else {
		cerr << "Unsupported mode " << mode << endl;
		exit(1);
	}

	this->tmp_closed = false;
	this->header_over = false;
	this->indexed = false;
	this->footer = nullptr;
	this->footer_discovery_ended = true;

	// Write the signature and the version at the beginning of the file
	if (this->is_writer) {
		uint8_t default_encoding = 0b00011110;
		// Signature
		uint8_t buff[] = {	'K', 'F', 'F',
							KFF_VERSION_MAJOR, KFF_VERSION_MINOR,
							default_encoding,
							0 /*uniqueness*/, 0 /*canonicity*/
						};

		this->write(buff, 8);

		this->indexed = true;
		this->end_position = 0;
	}
	// Read the header
	else if (this->is_reader) {
		// Header integrity marker
		uint8_t buff[4];
		this->read(buff, 3);
		if (buff[0] != 'K' or buff[1] != 'F' or buff[2] != 'F') {
			cerr << "Absent KFF signature at the beginning of the file." << endl;
			cerr << "Please check that the file is not corrupted" << endl;
			throw "Absent signature at the beginning";
		}

		// Version reading
		this->read(&this->major_version, 1);
		this->read(&this->minor_version, 1);
		if (KFF_VERSION_MAJOR < this->major_version or (KFF_VERSION_MAJOR == this->major_version and KFF_VERSION_MINOR < this->minor_version)) {
			cerr << "The software version " << (uint)KFF_VERSION_MAJOR << "." << (uint)KFF_VERSION_MINOR << " can't read files writen in version " << (uint)this->major_version << "." << (uint)this->minor_version << endl;
			throw "Unexpected version number";
		}
		// Encoding load
		this->read_encoding();
		// Read global flags
		uint8_t flag;
		this->read(&flag, 1);
		this->uniqueness = flag != 0;

		this->read(&flag, 1);
		this->canonicity = flag != 0;
		// Read metadata size
		this->read(buff, 4);
		load_big_endian(buff, 4, this->metadata_size);


		// Footer integrity marker
		unsigned long saved_position = this->tellp();
		this->jump_to(3, true);
		this->end_position = this->tellp();
		this->read(buff, 3);
		this->jump_to(saved_position);
		if (buff[0] != 'K' or buff[1] != 'F' or buff[2] != 'F') {
			cerr << "Absent KFF signature at the end of the file." << endl;
			cerr << "Please check that the file is not corrupted" << endl;
			throw "Absent signature at the end";
		}

		// Back to the start
		this->footer_discovery_ended = false;
		// Discover footer
		this->footer_discovery();
		this->index_discovery();
	}
}

void Kff_file::close(bool write_buffer) {
	if (this->is_writer) {
		// Write the index
		if (this->indexed)
			this->write_footer();
		// Write the signature
		char signature[] = {'K', 'F', 'F'};
		this->write((uint8_t *)signature, 3);
		
		// Write the end of the file
		if (write_buffer) {
			// The file was never opened
			if (not this->writing_started) {
				this->writing_started = true;
				this->fs.open(this->filename, fstream::binary | fstream::out);
			} else if (this->tmp_closed) {
				this->reopen();
			}
			// Write the buffer
			this->fs.write((char *)this->file_buffer, this->next_free);
			if (this->fs.fail()) {
				cerr << "Filesystem problem during buffer disk saving" << endl;
				exit(1);
			}
			this->file_size += this->next_free;
			this->next_free = 0;
		} else {
			this->delete_on_destruction = true;
		}

		// cout << "delete_on_destruction " << delete_on_destruction << endl;
		// cout << this->filename << endl;
		if (this->fs.is_open())
			this->fs.close();
	}
	else if (this->is_reader) {
		
	}

	this->tmp_closed = false;
	this->is_writer = false;
	this->is_reader = false;
}


Kff_file::~Kff_file() {
	this->close();

	delete[] this->file_buffer;
	if (this->delete_on_destruction and this->file_size > 0) {
		remove(this->filename.c_str());
	}

	if (this->footer != nullptr)
		delete this->footer;

	for (Section_Index * si : this->index)
		delete si;
}


void Kff_file::set_indexation(bool indexed) {
	if (this->is_writer)
		this->indexed = indexed;
}


void Kff_file::register_position(char section_type) {
	if (this->is_writer and this->indexed) {
		this->section_positions[this->tellp()] = section_type;
	}
}


void Kff_file::complete_header() {
	if (this->header_over)
		return;

	// If the metadata has not been read, jump over
	if (this->is_reader) {
		this->jump(this->metadata_size);
	}

	// If metadata has not been write, write a 0 byte one.
	else if (this->is_writer) {
		this->write_metadata(0, nullptr);
	}

	this->header_over = true;
}


void Kff_file::footer_discovery() {
	long current_pos = this->tellp();

	// Look at the footer
	this->jump_to(23, true);
	// Try to extract the footer size
	stringstream ss;
	char c = 'o';
	for (uint i=0 ; i<11 ; i++) {
		this->read((uint8_t *)&c, 1);
		ss << c;
	}
	if (ss.str().compare("footer_size") != 0) {
		return;
	}
	this->jump(1); // remove the '\0'

	uint64_t size = 0;
	uint8_t buff[8];
	this->read(buff, 8);
	load_big_endian(buff, 8, size);
	// Jump to value section start
	this->jump_to(size+3, true);

	this->footer = new Section_GV(this);
	this->footer->close();
	this->footer_discovery_ended = true;

	this->jump_to(current_pos);
}


void Kff_file::index_discovery() {
	long current_pos = this->tellp();
	bool header_over = this->header_over;
	this->complete_header();

	// Search in footer
	if (this->footer != nullptr and this->footer->vars.find("first_index") != this->footer->vars.end()) {
		this->indexed = true;
		this->read_index((long)this->footer->vars["first_index"]);
	}

	// Search first section
	if (not this->indexed) {
		char type = this->fs.peek();
		if (type == 'i') {
			this->indexed = true;
			this->read_index(this->tellp());
		}

	}

	this->header_over = header_over;
	this->index_discovery_ended = true;

	this->jump_to(current_pos);
}


void Kff_file::read_index(long position) {
	long init_pos = this->tellp();

	while (position != 0) {
		// Move to the beginning
		this->jump_to(position);
		// read the local index content
		Section_Index * si = new Section_Index(this);
		this->index.push_back(si);
		si->close();
		// Update index position to the next index section
		if (si->next_index == 0)
			position = 0;
		else {
			position = this->tellp() + si->next_index;
		}
	}

	this->jump_to(init_pos);
}


void Kff_file::read(uint8_t * bytes, unsigned long size) {
	if (not this->is_reader) {
		cerr << "Cannot read a file in writing mode." << endl;
		exit(1);
	}

	// Read in the file
	if (this->current_position < this->file_size) {
		// Read the end of the file and the beginning of the buffer
		if (this->current_position + size > this->file_size) {
			uint64_t fs_read_size = this->file_size - this->current_position;
			this->read(bytes, fs_read_size);
			this->read(bytes + fs_read_size, size - fs_read_size);
			return;
		}
		// Read inside of the file
		else {
			// File not opened
			if (not this->fs.is_open())
				this->fs.open(this->filename, fstream::binary | fstream::in);

			// long tp = this->fs.tellp();
			this->fs.read((char *)bytes, size);
			if (this->fs.fail()) {
				// cout << tp << endl;
				cerr << "Impossible to read the file " << this->filename << " on disk." << endl;
				exit(1);
			}
		}
	}
	// Read in the buffer
	else {
		// Compute the buffer positions to read
		uint64_t buffer_position = this->current_position - this->file_size;
		if (buffer_position + size > this->next_free) {
			cerr << "Read out of the file, Byte " << (this->file_size + this->next_free) << endl;
			exit(1);
		}

		memcpy(bytes, this->file_buffer + buffer_position, size);
	}
	
	this->current_position += size;
}

void Kff_file::write(const uint8_t * bytes, unsigned long size) {
	if (not this->is_writer) {
		if (this->is_reader)
			cerr << "Cannot write a file in reading mode." << endl;
		else
			cerr << "Cannot write a closed file" << endl;
		exit(1);
	}

	unsigned long buff_space = this->buffer_size - this->next_free;

	// Resize buffer
	while (buff_space < size and this->buffer_size < this->max_buffer_size) {
		// Enlarge the buffer
		this->buffer_size *= 2;
		uint8_t * next_buffer = new uint8_t[this->buffer_size];
		// Copy the previous values
		memcpy(next_buffer, this->file_buffer, this->next_free);
		buff_space = this->buffer_size - this->next_free;
		// Fill the empty part with 0s
		memset(next_buffer + this->next_free, 0, buff_space);
		// remove the previous space
		delete[] this->file_buffer;
		this->file_buffer = next_buffer;
	}

	// fill the buffer
	if (buff_space >= size) {
		memcpy(this->file_buffer + this->next_free, bytes, size);
		this->next_free += size;
	}
	// Not enought space, write the file
	else {
		// Open the file if needed
		if (not this->writing_started) {
			this->fs.open(this->filename, fstream::binary | fstream::out);
			this->writing_started = true;
		} else if (this->tmp_closed) {
			this->reopen();
		}

		this->fs.write((char*)this->file_buffer, this->next_free);
		this->fs.write((char*)bytes, size);
		this->file_size += this->next_free + size;
		this->next_free = 0;

		if (this->fs.fail()) {
			cerr << "File system error while writing " << this->filename << endl;
			exit(1);
		}
	}

	this->current_position += size;
}

void Kff_file::write_at(const uint8_t * bytes, unsigned long size, unsigned long position) {
	if (not this->is_writer) {
		if (this->is_reader)
			cerr << "Cannot write a file in reading mode." << endl;
		else
			cerr << "Cannot write a closed file" << endl;
		exit(1);
	}

	if (position > this->file_size + this->next_free) {
		cerr << "Cannot write after the last byte of the file." << endl;
		exit(1);
	}

	// Write the file on disk
	if (position < this->file_size) {
		// Only in file
		if (position + size <= this->file_size) {
			if (this->tmp_closed) {
				this->reopen();
			}
			this->fs.seekp(position);
			this->fs.write((char*)bytes, size);
			if (this->fs.fail()) {
				cerr << "File system error while writing " << this->filename << " at position " << position << endl;
				exit(1);
			}
			this->fs.seekp(this->file_size);
		}
		// On both file and buffer
		else {
			unsigned long in_file_size = this->file_size - position;
			// Write the file part
			this->write_at(bytes, in_file_size, position);
			// Write the buffer part
			this->write_at(bytes + in_file_size, size - in_file_size, position + in_file_size);
		}
	}
	// Write the buffer in RAM
	else {
		unsigned long corrected_position = position - this->file_size;
		
		// Write in the current buffer space
		if (corrected_position + size <= this->next_free) {
			memcpy(this->file_buffer + corrected_position, bytes, size);
		}
		// Spillover the buffer
		else {
			this->next_free = corrected_position;
			this->write(bytes, size);
		}
	}
}

unsigned long Kff_file::tellp() {
	return this->current_position;
}


void Kff_file::jump(long size) {
	// cout << "Jump " << this->current_position << " " << size << " / " << this->file_size << " " << this->next_free << endl;
	this->jump_to(this->current_position + size);
}

void Kff_file::jump_to(unsigned long position, bool from_end) {
	if (this->file_size + this->next_free < position) {
		cerr << "Jump out of the file." << endl;
		exit(1);
	}

	// Determine absolute position
	if (from_end) {
		position = this->file_size + this->next_free - position;
	}
	// cout << "position " << position << endl;

	// Jump into the written file
	if (position < this->file_size) {
		this->fs.seekp(position);
	}
	// Jump into the buffer
	else /*if (this->current_position < this->file_size)*/ {
		this->fs.seekg(0, this->fs.end);
	}
	this->current_position = position;
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


void Kff_file::write_footer() {
	Section_Index si(this);

	// Compute end position
	long position = si.beginning + 17 + 9 * this->section_positions.size();
	// Add the values
	for (auto & it : this->section_positions) {
		si.register_section(it.second, it.first - position);
	}

	si.close();

	// Write a value section to register everything
	Section_GV sgv(this);
	sgv.write_var("first_index", si.beginning);
	sgv.write_var("footer_size", 9 + 2 * (12 + 8));
	sgv.close();
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
	this->write_at(&code, 1, 5);
}

void Kff_file::set_uniqueness(bool uniqueness) {
	uint8_t bit_uniq = uniqueness ? 1 : 0;
	this->write_at(&bit_uniq, 1, 6);
}
void Kff_file::set_canonicity(bool canonicity) {
	uint8_t bit_canon = canonicity ? 1 : 0;
	this->write_at(&bit_canon, 1, 7);
}

void Kff_file::write_encoding(uint8_t * encoding) {
	this->write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
}

void Kff_file::read_encoding() {
	uint8_t code, a, c, g, t;
	// Get code values
	this->read(&code, 1);

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
	if (this->header_over) {
		cerr << "The metadata have to be written prior to other content." << endl;
		exit(1);
	}

	uint8_t buff[4];
	store_big_endian(buff, 4, size);
	this->write(buff, 4);
	this->write(data, size);

	this->header_over = true;
}


void Kff_file::read_metadata(uint8_t * data) {
	this->read(data, this->metadata_size);
	this->header_over = true;
}

bool Kff_file::jump_next_section() {
	if (not is_reader)
		return false;
	char section_type = read_section_type();
	if (this->current_position == this->file_size + this->next_free)
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

	if (this->current_position < this->file_size) {
		return this->fs.peek();
	}
	else {
		return (char)this->file_buffer[this->current_position - this->file_size];
	}
}


Section::Section(Kff_file * file) {
	this->file = file;

	if (not file->header_over and file->footer_discovery_ended) {
		file->complete_header();
	}

	this->beginning = file->tellp();
}

void Section::close() {
	this->file = nullptr;
}


Section * SectionBuilder::build(Kff_file * file) {
	char type = file->read_section_type();
	switch (type) {
		case 'i':
			return new Section_Index(file);
		case 'v':
			return new Section_GV(file);
		case 'r':
			return new Section_Raw(file);
		case 'm':
			return new Section_Minimizer(file);
		default:
			cerr << "Unknown section " << type << "(" << (uint)type << ")" << endl;
			throw std::runtime_error("Unknown section type" + type);
	}
}


// ----- Global variables sections -----

Section_GV::Section_GV(Kff_file * file) : Section(file) {
	this->nb_vars = 0;
	this->file->global_vars.clear();

	if (this->file->is_reader) {
		this->read_section();
	}

	if (file->is_writer) {
		if (file->indexed)
			file->register_position('v');
		char type = 'v';
		this->file->write((uint8_t *)&type, 1);
	}
}

void Section_GV::write_var(const string & var_name, uint64_t value) {
	this->nb_vars += 1;
	this->vars[var_name] = value;
	this->file->global_vars[var_name] = value;
}

void Section_GV::read_section() {
	char type = '\0';
	this->file->read((uint8_t *)&type, 1);
	if (type != 'v')
		throw "The section do not start with the 'v' char, you can't open a Global Variable section.";

	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_vars);
	for (uint64_t i=0 ; i<nb_vars ; i++) {
		this->read_var();
	}
}

void Section_GV::read_var() {
	if (file->tellp() >= file->end_position)
		throw "eof reached before the end of the variable section";

	// Name reading
	stringstream ss;
	char c = 'o';
	this->file->read((uint8_t *)&c, 1);
	while (c != '\0') {
		ss << c;
		this->file->read((uint8_t *)&c, 1);
	}

	// Value reading
	uint64_t value = 0;
	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, value);

	// Saving
	string name = ss.str();
	this->vars[name] = value;
	this->file->global_vars[name] = value;
}

void Section_GV::copy(Kff_file * file) {
	// Remove empty variable sections
	if (this->vars.size() == 0)
		return;

	// Open the copy
	Section_GV sgv(file);
	// Copy all the variables
	for (const auto & it : this->vars) {
		sgv.write_var(it.first, it.second);
	}
	// Clos the copy
	sgv.close();
}

void Section_GV::close() {
	if (file->is_writer) {
		uint8_t buff[8];
		// write the number of block values
		store_big_endian(buff, 8, this->nb_vars);
		this->file->write(buff, 8);
		// Write the variables
		for (std::map<std::string,uint64_t>::iterator var_tuple=this->vars.begin() ; var_tuple != this->vars.end() ; var_tuple++) {
			const string & name = var_tuple->first;
			this->file->write((const uint8_t *)name.c_str(), name.length()+1);
			store_big_endian(buff, 8, var_tuple->second);
			this->file->write(buff, 8);
		}
	}

	Section::close();
}



Section_Index::Section_Index(Kff_file * file) : Section(file) {
	char type;
	uint8_t buff[8];

	this->next_index = 0;

	if (this->file->is_reader) {
		this->file->read((uint8_t *)&type, 1);
		if (type != 'i')
			throw "The section do not start with the 'i' char, you can not open an Index section.";

		uint64_t nb_vars;
		this->file->read(buff, 8);
		load_big_endian(buff, 8, nb_vars);
		for (uint64_t i=0 ; i<nb_vars ; i++) {
			int64_t idx = 0;
			this->file->read((uint8_t *)&type, 1);
			this->file->read(buff, 8);
			load_big_endian(buff, 8, idx);
			this->index[idx] = type;
		}

		if (nb_vars != this->index.size())
			throw "index collision in i section";

		this->file->read(buff, 8);
		load_big_endian(buff, 8, this->next_index);
	}
}

void Section_Index::register_section(char section_type, int64_t pos) {
	this->index[pos] = section_type;
}

void Section_Index::set_next_index(int64_t index) {
	this->next_index = index;
}

// void Section_Index::copy(Kff_file * file) {
// 	cerr << "You are trying to copy an index from a file to another. ";
// 	cerr << "As the positions can be different between the files, this operation is not allowed." << endl;
// 	exit(2);
// }

void Section_Index::close() {
	if (this->file->is_writer) {
		uint8_t buff[8];
		// Section header
		char type = 'i';
		this->file->write((uint8_t *)&type, 1);
		store_big_endian(buff, 8, this->index.size());
		this->file->write(buff, 8);
		// Write index
		for (std::map<int64_t, char>::iterator it=this->index.begin(); it!=this->index.end(); ++it) {
			// Section type
			type = it->second;
			this->file->write((uint8_t *)&type, 1);
			// Section index
			store_big_endian(buff, 8, it->first);
			this->file->write(buff, 8);
		}
		store_big_endian(buff, 8, this->next_index);
		this->file->write(buff, 8);
	}

	Section::close();
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

Section_Raw::Section_Raw(Kff_file * file) : Section(file){
	if (file->global_vars.find("k") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing k variable";
	if(file->global_vars.find("max") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing max variable";
	if(file->global_vars.find("data_size") == file->global_vars.end())
		throw "Impossible to read the raw section due to missing data_size variable";
	
	uint64_t k = file->global_vars["k"];
	uint64_t max = file->global_vars["max"];
	uint64_t data_size = file->global_vars["data_size"];

	this->nb_blocks = 0;

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
		if (file->indexed)
			file->register_position('r');
		char type = 'r';
		this->file->write((uint8_t *)&type, 1);
		this->file->write((uint8_t *)&this->nb_blocks, 8);
	}
}

uint32_t Section_Raw::read_section_header() {
	char type;
	this->file->read((uint8_t *)&type, 1);
	if (type != 'r')
		throw "The section do not start with the 'r' char, you can't open a Raw sequence section.";

	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_blocks);
	this->remaining_blocks = this->nb_blocks;

	return this->nb_blocks;
}

void Section_Raw::write_compacted_sequence(uint8_t* seq, uint64_t seq_size, uint8_t * data_array) {
	uint8_t buff[8];
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size - k + 1;
	store_big_endian(buff, this->nb_kmers_bytes, nb_kmers);
	this->file->write(buff, this->nb_kmers_bytes);
	// 2 - Write sequence
	uint64_t seq_bytes_needed = (seq_size + 3) / 4;
	this->file->write(seq, seq_bytes_needed);
	// 3 - Write data
	uint64_t data_bytes_needed = data_size * nb_kmers;
	this->file->write(data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Raw::read_compacted_sequence(uint8_t* seq, uint8_t* data) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	this->file->read(seq, seq_bytes_needed);
	// 3 - Read the data.
	uint64_t data_bytes_used = data_size * nb_kmers_in_block;
	this->file->read(data, data_bytes_used);

	this->remaining_blocks -= 1;

	return nb_kmers_in_block;
}


void Section_Raw::copy(Kff_file * file) {
	uint max_nucl = this->k + this->max - 1;
	uint8_t * seq_buffer = new uint8_t[(max_nucl + 3) / 4];
	uint8_t * data_buffer = new uint8_t[this->max * this->data_size];

	// Open the copy
	Section_Raw sr(file);
	// Copy all the variables
	for (uint i=0 ; i<this->nb_blocks ; i++) {
		// Read
		uint64_t size = this->read_compacted_sequence(seq_buffer, data_buffer);
		// Rewrite
		sr.write_compacted_sequence(seq_buffer, this->k + size - 1, data_buffer);
	}
	// Clos the copy
	sr.close();
	delete[] seq_buffer;
	delete[] data_buffer;
}


void Section_Raw::jump_sequence() {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	// 3 - Determine the data size
	size_t data_bytes_used = data_size * nb_kmers_in_block;
	// 4 - Jumb over the 
	file->jump(seq_bytes_needed + data_bytes_used);
	this->remaining_blocks -= 1;
}


void Section_Raw::close() {
	if (this->file->is_writer) {
		uint8_t buff[8];
		store_big_endian(buff, 8, this->nb_blocks);
		this->file->write_at(buff, 8, this->beginning + 1);
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}

	Section::close();
}



// ----- Minimizer sequence section -----

Section_Minimizer::Section_Minimizer(Kff_file * file) : Section(file) {
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

	this->nb_blocks = 0;
	this->remaining_blocks = 0;

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
		if (file->indexed)
			file->register_position('m');

		char type = 'm';
		this->file->write((uint8_t *)&type, 1);
		this->file->write(this->minimizer, this->nb_bytes_mini);
		this->file->write((uint8_t *)&this->nb_blocks, 8);
	}
}

Section_Minimizer::~Section_Minimizer() {
	delete[] this->minimizer;
}

Section_Minimizer& Section_Minimizer::operator= ( Section_Minimizer && sm) {
	file = sm.file;
	sm.file = nullptr;
	beginning = sm.beginning;
	nb_blocks = sm.nb_blocks;

	m = sm.m;
	k = sm.k;
	max = sm.max;
	data_size = sm.data_size;

	this->remaining_blocks = sm.remaining_blocks;
	this->nb_kmers_bytes = nb_kmers_bytes;

	nb_bytes_mini = sm.nb_bytes_mini;
	std::swap(minimizer, sm.minimizer);

	return *this;
}

void Section_Minimizer::write_minimizer(uint8_t * minimizer) {
	this->file->write_at(minimizer, this->nb_bytes_mini, this->beginning+1);

	// uint64_t pos = file->fs.tellp();
	// file->fs.seekp(this->beginning+1);
	// file->fs.write((char *)minimizer, this->nb_bytes_mini);
	// memcpy(this->minimizer, minimizer, this->nb_bytes_mini);
	// file->fs.seekp(pos);
}

uint32_t Section_Minimizer::read_section_header() {
	// Verify section type
	char type;
	this->file->read((uint8_t *)&type, 1);
	if (type != 'm')
		throw "The section do not start with the 'm' char, you can't open a Minimizer sequence section.";

	// Read the minimizer
	this->file->read(this->minimizer, this->nb_bytes_mini);

	// Read the number of following blocks
	uint8_t buff[8];
	this->file->read(buff, 8);
	load_big_endian(buff, 8, this->nb_blocks);
	this->remaining_blocks = this->nb_blocks;
	return nb_blocks;
}

void Section_Minimizer::write_compacted_sequence_without_mini(uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array) {
	uint8_t buff[8];
	// 1 - Write nb kmers
	uint64_t nb_kmers = seq_size + m - k + 1;
	store_big_endian(buff, this->nb_kmers_bytes, nb_kmers);
	file->write(buff, this->nb_kmers_bytes);
	// 2 - Write minimizer position
	store_big_endian(buff, this->mini_pos_bytes, mini_pos);
	file->write(buff, this->mini_pos_bytes);
	// 3 - Write sequence with chopped minimizer
	uint64_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	this->file->write(seq, seq_bytes_needed);
	// 4 - Write data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers);
	this->file->write(data_array, data_bytes_needed);

	this->nb_blocks += 1;
}

uint64_t Section_Minimizer::read_compacted_sequence_without_mini(uint8_t* seq, uint8_t* data, uint64_t & mini_pos) {
	uint8_t buff[8];
	uint64_t nb_kmers_in_block = 1;
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 2 - Read the minimizer position
	file->read(buff, this->mini_pos_bytes);
	load_big_endian(buff, this->mini_pos_bytes, mini_pos);
	// 3 - Read the sequence
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = bytes_from_bit_array(2, seq_size);
	file->read(seq, seq_bytes_needed);
	// 4 - Read the data
	uint64_t data_bytes_needed = bytes_from_bit_array(data_size*8, nb_kmers_in_block);
	file->read(data, data_bytes_needed);
	// cout << data_bytes_needed << endl;

	this->remaining_blocks -= 1;
	return nb_kmers_in_block;
}

void Section_Minimizer::copy(Kff_file * file) {
	uint max_nucl = this->k + this->max - 1;
	uint8_t * seq_buffer = new uint8_t[(max_nucl + 3) / 4];
	uint8_t * data_buffer = new uint8_t[this->max * this->data_size];
	uint64_t mini_pos = 0;

	// Open the copy
	Section_Minimizer sm(file);
	sm.write_minimizer(this->minimizer);

	// Copy all the variables
	for (uint i=0 ; i<this->nb_blocks ; i++) {
		// Read
		uint64_t size = this->read_compacted_sequence_without_mini(seq_buffer, data_buffer, mini_pos);
		// Rewrite
		sm.write_compacted_sequence_without_mini(seq_buffer, this->k + size - 1 - this->m, mini_pos, data_buffer);
	}
	// Clos the copy
	sm.close();
	delete[] seq_buffer;
	delete[] data_buffer;
}

void Section_Minimizer::jump_sequence() {
	uint64_t nb_kmers_in_block = 1;
	uint8_t buff[8];
	// 1 - Read the number of kmers in the sequence
	if (nb_kmers_bytes != 0) {
		file->read(buff, this->nb_kmers_bytes);
		// 2 - Convert from big endian
		load_big_endian(buff, this->nb_kmers_bytes, nb_kmers_in_block);
	}
	// 3 - Determine the sequence size
	size_t seq_size = nb_kmers_in_block + k - m - 1;
	size_t seq_bytes_needed = (seq_size + 3) / 4;
	// 3 - Determine the data size
	size_t data_bytes_used = data_size * nb_kmers_in_block;
	// 4 - Jumb over the 
	file->jump((long)(this->mini_pos_bytes + seq_bytes_needed + data_bytes_used));

	this->remaining_blocks -= 1;
}


/* Bitshift to the left all the bits in the array with a maximum of 7 bits.
 * Overflow on the left will be set into the previous cell.
 */
static void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	if (length > 0) {
		for (uint64_t i=0 ; i<length-1 ; i++) {
			bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
		}
		bitarray[length-1] <<= bitshift;
	}
}

/* Similar to the previous function but on the right */
static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	if (length > 0) {
		for (uint64_t i=length-1 ; i>0 ; i--) {
			bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
		}
		bitarray[0] >>= bitshift;
	}
}

/* Fusion to bytes into one.
 * The merge_index higher bits are from left_bits the others from right_bits
 */
static uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}

void Section_Minimizer::write_compacted_sequence (uint8_t* seq, uint64_t seq_size, uint64_t mini_pos, uint8_t * data_array) {
	// Compute all the bit and byte quantities needed.
	uint64_t seq_bytes = bytes_from_bit_array(2, seq_size);
	uint left_offset_nucl = (4 - seq_size % 4) % 4;

	// 1 - Prepare space for sequence manipulation
	uint8_t * seq_copy = new uint8_t[seq_bytes];
	memcpy(seq_copy, seq, seq_bytes);
	
	// 2 - Move the suffix to the bytes where the minimiser started
	uint mini_start_byte = (mini_pos + left_offset_nucl) / 4;
	uint suff_start_byte = (mini_pos + m + left_offset_nucl) / 4;
	uint suff_bytes = seq_bytes-suff_start_byte;
	for (uint i=0 ; i<suff_bytes ; i++) {
		seq_copy[mini_start_byte + i] = seq_copy[suff_start_byte + i];
	}

	// 3 - shift the suffix to align with the minimizer position
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
	uint64_t seq_left_offset = (4 - (seq_size % 4)) % 4;
	uint64_t no_mini_size = seq_size - m;
	uint64_t no_mini_bytes = bytes_from_bit_array(2, no_mini_size);
	uint64_t no_mini_left_offset = (4 - (no_mini_size % 4)) % 4;
	// Shift the whole sequence to the left to have np padding on byte 0.
	leftshift8(seq, no_mini_bytes, no_mini_left_offset*2);


	// Prepare the suffix
	uint8_t * suffix = new uint8_t[seq_bytes];
	memset(suffix, 0, seq_bytes);
	uint suff_nucl = seq_size - m - mini_pos;
	// Values inside of seq before any change
	uint no_mini_suff_start_nucl = mini_pos;
	uint no_mini_suff_start_byte = no_mini_suff_start_nucl / 4;
	uint no_mini_suff_bytes = no_mini_bytes - no_mini_suff_start_byte;
	memcpy(suffix, seq + no_mini_suff_start_byte, no_mini_suff_bytes);
	// Shift to the left
	uint no_mini_suff_offset = no_mini_suff_start_nucl % 4;
	leftshift8(suffix, no_mini_suff_bytes, no_mini_suff_offset * 2);



	// Prepare the minimizer
	uint8_t * mini = new uint8_t[seq_bytes];
	memset(mini, 0, seq_bytes);
	memcpy(mini, this->minimizer, nb_bytes_mini);
	// Shift to the left
	uint mini_offset = (4 - (m % 4)) % 4;
	leftshift8(mini, nb_bytes_mini, mini_offset * 2);


	// Align the minimizer
	uint final_mini_start_nucl = mini_pos;
	uint final_mini_start_byte = mini_pos / 4;
	uint final_mini_offset = final_mini_start_nucl % 4;
	uint final_mini_byte_size = (m + final_mini_offset + 3) / 4;
	rightshift8(mini, seq_bytes, final_mini_offset * 2);

	// Align the suffix with the end of the minimizer
	uint final_suff_start_nucl = final_mini_start_nucl + m;
	uint final_suff_start_byte = final_suff_start_nucl / 4;
	uint final_suff_offset = final_suff_start_nucl % 4;
	uint final_suff_byte_size = (suff_nucl + final_suff_offset + 3) / 4;
	rightshift8(suffix, seq_bytes, final_suff_offset * 2);


	// Merge minimizer
	seq[final_mini_start_byte] = fusion8(
		seq[final_mini_start_byte],
		mini[0],
		final_mini_offset * 2
	);
	for (uint idx=1 ; idx<final_mini_byte_size ; idx++) {
		seq[final_mini_start_byte+idx] = mini[idx];
	}

	// Merge the suffix
	seq[final_suff_start_byte] = fusion8(
		seq[final_suff_start_byte],
		suffix[0],
		final_suff_offset * 2
	);
	for (uint64_t idx=1 ; idx<final_suff_byte_size ; idx++) {
		seq[final_suff_start_byte+idx] = suffix[idx];
	}

	// Align everything to the right
	rightshift8(seq, seq_bytes, seq_left_offset * 2);

	delete[] suffix;
	delete[] mini;
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
		uint8_t tmp[8];
		store_big_endian(tmp, 8, this->nb_blocks);
		this->file->write_at(tmp, 8, this->beginning + 1l + (long)this->nb_bytes_mini);
		// // Save current position
		// long position = this->file.tellp();
		// // Go write the number of variables in the correct place
		// fs.seekp(this->beginning + 1l + (long)this->nb_bytes_mini);
		// write_value(nb_blocks, fs);
		// fs.seekp(position);
	}

	if (file->is_reader) {
		// Jump over remaining sequences of the section
		while (this->remaining_blocks > 0)
			this->jump_sequence();
	}

	Section::close();
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

	this->k = 0;
	this->max = 0;
	this->data_size = 0;

	this->has_next();
}

Kff_reader::~Kff_reader() {
	delete[] this->current_kmer;
	delete[] this->current_data;
	if (this->current_section != NULL)
		delete this->current_section;

	for (uint i=0 ; i<4 ; i++)
		delete[] this->current_shifts[i];
	delete[] this->current_shifts;

	delete this->file;
}

void Kff_reader::read_until_first_section_block() {
	while (current_section == NULL or remaining_blocks == 0) {
		if (this->file->tellp() == this->file->end_position) {
			break;
		}

		// char section_type = this->file->read_section_type();
		char section_type = file->read_section_type();
		// --- Update data structure sizes ---
		if (section_type == 'v') {
			// Read the global variable block
			Section_GV gvs(file);
			// Update sequence size if k or max change
			if (gvs.vars.find("k") != gvs.vars.end()
				or gvs.vars.find("max") != gvs.vars.end()) {
				// Compute the max size of a sequence
				this->k = this->file->global_vars["k"];
				this->max = this->file->global_vars["max"];
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
				this->data_size = data_size;
				this->max = this->file->global_vars["max"];
				uint64_t max_size = data_size * max;
				delete[] this->current_data;
				this->current_data = new uint8_t[max_size];
				memset(this->current_data, 0, max_size);
			}
		}
		// Mount data from the files to the datastructures.
		else if (section_type == 'i') {
			Section_Index index(file);
			index.close();
		} else {
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
	for (uint8_t i=1 ; i<min((uint64_t)4, remaining_kmers) ; i++) {
		// Copy
		memcpy(current_shifts[i], current_sequence, current_seq_bytes);
		// Shift
		rightshift8(current_shifts[i], current_seq_bytes, 2 * i);
	}
}

bool Kff_reader::has_next() {
	if (current_section == NULL and (file->end_position > file->tellp()))
		read_until_first_section_block();
	return file->end_position > file->tellp();
}

uint64_t Kff_reader::next_block(uint8_t* & sequence, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		sequence = NULL;
		data = NULL;
		return 0;
	}

	uint64_t nb_kmers = current_section->read_compacted_sequence(sequence, data);
	
	remaining_kmers = 0;
	remaining_blocks -= 1;
	if (remaining_blocks == 0) {
		delete current_section;
		current_section = NULL;
	}

	return nb_kmers;
}

bool Kff_reader::next_kmer(uint8_t* & kmer, uint8_t* & data) {
	// Verify the abylity to find another kmer in the file.
	if (!this->has_next()){
		kmer = NULL;
		data = NULL;
		return false;
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
	uint64_t end_nucl = start_nucl + this->k - 1;
	uint64_t end_byte = end_nucl / 4;

	memcpy(current_kmer, current_shifts[right_shift]+start_byte, end_byte-start_byte+1);
	kmer = current_kmer;
	data = current_data + (current_seq_kmers - remaining_kmers) * this->data_size;
	
	// Read the next block if needed.
	remaining_kmers -= 1;
	if (remaining_kmers == 0) {
		remaining_blocks -= 1;
		if (remaining_blocks == 0) {
			delete current_section;
			current_section = NULL;
		}
	}

	return true;
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












