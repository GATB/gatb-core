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

    Range<u_int32_t>::Iterator itRange (1, nbIter);

    SubjectIterator<u_int32_t> it (itRange, nbIter/100);

    Progress progress (nbIter, "Allocating blocks");
    it.addObserver (progress);

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
    cout << "totalPhysMem = " << System::info().getMemoryPhysicalTotal() << endl;
    cout << "physMemUsed  = " << System::info().getMemoryPhysicalUsed()  << endl;
    cout << "buffersMem   = " << System::info().getMemoryBuffers()       << endl;

    /** Provided in MBytes */
    u_int64_t toErase = 1024 * 1024 * (argc >= 2 ? atol (argv[1]) : 0);

    clear (toErase);

    return EXIT_SUCCESS;
}
