#ifndef KFF_IO_MERGE
#define KFF_IO_MERGE


#include <vector>
#include <string>

#include "kff_io.hpp"


using namespace std;



void kff_merge2(const vector<Kff_file *> & files, string output);
void kff_merge (const vector<string> inputs, string output);

#endif
