#include <gatb/system/impl/System.hpp>
#include <ctime> // for time
#include <iostream> // for time (and maybe other things?)
#include <iomanip> // for cout mods


#include "logging.hpp"

namespace gatb { namespace core { namespace debruijn { namespace impl  {
bool bcalm_logging = true;

unsigned long logging(std::string message="")
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    if (bcalm_logging) 
    {
        std::cout << std::setiosflags(std::ios::right);
        std::cout << std::resetiosflags(std::ios::left);
        std::cout << std::setw(40) << std::left << message << "      ";
    }
    char tmp[128];
    snprintf (tmp, sizeof(tmp), "  %02d:%02d:%02d  ",
            now->tm_hour, now->tm_min, now->tm_sec);
    if (bcalm_logging) 
        std::cout << tmp ;

    // using Progress.cpp of gatb-core
    u_int64_t mem = gatb::core::system::impl::System::info().getMemorySelfUsed() / 1024;
    u_int64_t memMaxProcess = gatb::core::system::impl::System::info().getMemorySelfMaxUsed() / 1024;
    snprintf (tmp, sizeof(tmp), "   memory [current, maxRSS]: [%4lu, %4lu] MB ",
            mem, memMaxProcess);

    if (bcalm_logging) 
        std::cout << tmp << std::endl;
    return mem;
}


}}}}
