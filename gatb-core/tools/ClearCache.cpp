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
#include <gatb/tools/misc/api/Range.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <iostream>
#include <vector>

using namespace std;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace gatb::core::tools::dp::impl;

/********************************************************************************/
void clear (u_int64_t toErase)
{
    u_int64_t totalPhysical = System::info().getMemoryPhysicalTotal();

    if (toErase == 0)  {  toErase = totalPhysical; }

    size_t blockSize = 32*1024;

    u_int32_t nbIter = toErase/blockSize;

    vector<void*> buffers (nbIter, 0);

    SubjectIterator<u_int32_t> it (new Range<u_int32_t>::Iterator (1, nbIter), nbIter/100);

    it.addObserver (new ProgressTimer (nbIter, "clear cache"));

    size_t i=0;
    for (it.first(); !it.isDone(); it.next())
    {
        void* b = calloc (blockSize, 1);
        if (b != 0)  {  buffers[i]=b; i++; }
    }

    //#if 0
    //    size_t j=0;
    //    for (size_t k=0; k<nbIter; k++)
    //    {
    //        void* b = buffers[k];
    //        if (b != 0)  { free (b); j++; }
    //    }
    //#endif
}

/********************************************************************************/

int main (int argc, char* argv[])
{
    //cout << "totalPhysMem = " << System::info().getMemoryPhysicalTotal() << endl;
    //cout << "physMemUsed  = " << System::info().getMemoryPhysicalUsed()  << endl;
    //cout << "buffersMem   = " << System::info().getMemoryBuffers()       << endl;

    /** Provided in MBytes */
    u_int64_t toErase = 1024 * 1024 * (argc >= 2 ? atol (argv[1]) : 0);

    clear (toErase);

    return EXIT_SUCCESS;
}
