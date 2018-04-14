/* trying to write a large bagfile
	mimics bglue's prepare_uf
	attempts to reproduce a bug from seb raguideau, very large assembly
 * */

#include <chrono>
#define get_wtime() chrono::system_clock::now()
#define diff_wtime(x,y) chrono::duration_cast<chrono::nanoseconds>(y - x).count()


#include <gatb/system/impl/System.hpp>

#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/bank/impl/BankStrings.hpp>
#include <gatb/bank/impl/BankSplitter.hpp>
#include <gatb/bank/impl/BankRandom.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/bank/api/IBank.hpp>

#include <gatb/tools/collections/api/Bag.hpp>

#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>


#include <iostream>
#include <memory>

using namespace std;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::math;
using namespace gatb::core::tools::dp;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

typedef uint32_t partition_t;

void doit()
{

    BagFile<partition_t> * bagf = new BagFile<partition_t>( "temp.largefile"); LOCAL(bagf); 
    Bag<partition_t> * currentbag =  new BagCache<partition_t> (  bagf, 10000 ); LOCAL(currentbag);

    //uint64_t bound = 4600000000; // just over 32 bits
    uint64_t bound = 1389239905; // just over 4.2 GB file, what's being written by bglue's problematic instance
    for (uint64_t i = 0; i < bound; i++)
    {
        if (i % 100000000 == 0)
            std::cout << " wrote " << i << " elements" << std::endl;
	    partition_t cur = i;
	    currentbag->insert(cur);
    }
    
    std::cout << " flushing" << std::endl;
    currentbag->flush();
    std::cout << " done" << std::endl;
          
}
int main (int argc, char* argv[])
{
    try
    {
    	doit();
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }



}
