
#include <vector>
#include <string>
#include <iostream>

#include "kff_io.hpp"
#include "merge.hpp"


using namespace std;


void kff_merge(const vector<string> inputs, string output) {
	vector<Kff_file*> files;
	files.reserve(inputs.size());
	for(auto & input : inputs) {
    files.push_back(new Kff_file(input, "r"));  
	}

	kff_merge2(files, output);

    for (Kff_file* file : files)
    {
        delete file;
    }
}

void kff_merge2(const vector<Kff_file*> &files, string output) {
	// Useful variables
	const long buffer_size = 1048576; // 1 MB
	uint8_t buffer[1048576];
	uint8_t global_encoding[4];

	// Read the encoding of the first file and push it as outcoding
	for (uint i=0 ; i<4 ; i++)
		global_encoding[i] = files[0]->encoding[i];
	
	// Write header of the output
	Kff_file outfile(output, "w");
	outfile.set_indexation(true);
	outfile.set_uniqueness(true);
	outfile.set_canonicity(true);
	outfile.write_encoding(
		global_encoding[0],
		global_encoding[1],
		global_encoding[2],
		global_encoding[3]
	);
	// Set metadata
	std::string meta = "Merged file";
	outfile.write_metadata(meta.length(), (uint8_t *)meta.c_str());

	// remember index previous position for chaining
	long last_index = 0;
	// Footers
	map<string, uint64_t> footer_values;

	// Append each file one by one
	for (Kff_file * infile : files) {
		// Encoding verification
		for (uint i=0 ; i<4 ; i++) {
			if (infile->encoding[i] != global_encoding[i]) {
				cerr << "Wrong encoding for file " << infile->filename << endl;
				cerr << "Its nucleotide encoding is different from previous kff files." << endl;
				cerr << "Please first use 'kff-tools translate' to have the same encoding" << endl;
				exit(1);
			}
		}

		// NB: Automatic jump over metadata due to API
		// Read section by section
		char section_type = infile->read_section_type();
		while(infile->tellp() != infile->end_position) {
			vector<string> to_copy;
			long size, end_byte, begin_byte;

			switch (section_type) {
				// Write the variables that change from previous sections (possibly sections from other input files)
				case 'v':
				{
					unordered_map<string, uint64_t> variables;

					// Read variables
					while (section_type == 'v') {
						Section_GV sgv(infile);

						// Is it a footer ?
						if (sgv.vars.find("footer_size") != sgv.vars.end()) {
							for (auto& tuple : sgv.vars)
								if (tuple.first != "footer_size" and tuple.first != "first_index") {
									if (footer_values.find(tuple.first) == footer_values.end())
										footer_values[tuple.first] = tuple.second;
									else
										// Sum up the common footer values
										footer_values[tuple.first] += tuple.second;
								}
							break;
						}
						// Not a footer
						for (auto & p : sgv.vars)
							variables[p.first] = p.second;
						sgv.close();

						// Update section_type
						section_type = infile->read_section_type();
					}

					// Does it need a rewrite ?
					bool v_section_needed = false;
					for (auto & p : variables) {
						if (outfile.global_vars.find(p.first) == outfile.global_vars.end()
								or outfile.global_vars[p.first] != p.second) {
							v_section_needed = true;
							break;
						}
					}

					// cout << "V needed ?"
					// Rewrite
					if (v_section_needed) {
						Section_GV sgv(&outfile);
						for (auto & p : variables) {
							sgv.write_var(p.first, p.second);
						}
						sgv.close();
					}
				}
				break;

				// copy the sequence section from input to output
				case 'r':
				case 'm':
				// Analyse the section size
				begin_byte = infile->tellp();
				infile->jump_next_section();
				end_byte = infile->tellp();
				size = end_byte - begin_byte;
				infile->jump(-size);

				// Read from input and write into output
				while (size > 0) {
					size_t size_to_copy = size > buffer_size ? buffer_size : size;

					infile->read(buffer, size_to_copy);
					outfile.write(buffer, size_to_copy);

					size -= size_to_copy;
				}
				break;
				case 'i': {
				// read section and compute its size
				Section_Index si(infile);
				si.close();
				// long file_size = infile->tellp() - si.beginning - 8l;
				// infile->jump(-file_size - 8);

				// // Save the position in the file for later chaining
				// long i_position = outfile.tellp();
				// size = file_size;
				// // Copy section (except the chaining part)
				// // Read from input and write into output
				// while (size > 0) {
				// 	size_t size_to_copy = size > buffer_size ? buffer_size : size;

				// 	infile->read(buffer, size_to_copy);
				// 	outfile.write(buffer, size_to_copy);

				// 	size -= size_to_copy;
				// }
				// // Jump over the last value of infile
				// infile->jump(8);
				// // Chain the section and save its position
				// long i_relative = last_index - (i_position + file_size + 8l);
				// if (last_index == 0)
				// 	i_relative = 0;
				// for (uint i=0 ; i<8 ; i++) {
				// 	uint8_t val = (uint8_t)(i_relative >> (56 - 8 * i));
				// 	outfile.write(&val, 1);
				// }
				// // write_value(last_index, outfile.fs);
				// last_index = i_position;
				} break;

				default:
					cerr << "Unknown section type " << section_type << " in file " << infile->filename << endl;
					exit(2);
			}

			// Prepare next section
			section_type = infile->read_section_type();
		}

		infile->close();
	}

	// Write footer
	if (last_index != 0) {
		Section_GV sgv(&outfile);
		long size = 9;

		for (const auto & tuple : footer_values) {
			sgv.write_var(tuple.first, tuple.second);
			size += tuple.first.length() + 1 + 8;
		}
		sgv.write_var("first_index", last_index);
		sgv.write_var("footer_size", size + 2 * (12 + 8));
		sgv.close();

	}

	outfile.close();
}

