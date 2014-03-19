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

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/kmer/impl/Model.hpp>

#include <iostream>
#include <vector>

using namespace std;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::misc;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/

int main (int argc, char* argv[])
{

    if (argc < 2)
    {
        cerr << "you must provide at least 1 arguments. Arguments are:" << endl;
        cerr << "   1) uri of the file"  << endl;
        cerr << "   2) dump (0 default or 1)"  << endl;
        return EXIT_FAILURE;
    }

    string filename = argv[1];

    bool dump = (argc >= 3);

    u_int64_t nbItems  = 0;
    kmer_type checksum = 0;

    IteratorFile<kmer_type> it (filename);
    for (it.first(); !it.isDone(); it.next())
    {
        if (dump)  { cout << "0x" << *it << endl; }
        checksum = checksum + *it;
        nbItems ++;
    }

    cout << "FOUND " << nbItems << " WITH CHECKSUM " << checksum << endl;

    return EXIT_SUCCESS;
}
