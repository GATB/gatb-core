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

#include <gatb/system/impl/SystemInfoCommon.hpp>

#include <string.h>
#include <unistd.h>

using namespace std;

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

/*********************************************************************
                #        ###  #     #  #     #  #     #
                #         #   ##    #  #     #   #   #
                #         #   # #   #  #     #    # #
                #         #   #  #  #  #     #     #
                #         #   #   # #  #     #    # #
                #         #   #    ##  #     #   #   #
                #######  ###  #     #   #####   #     #
*********************************************************************/

#ifdef __linux__

#include "sys/sysinfo.h"

/********************************************************************************/
size_t SystemInfoLinux::getNbCores () const
{
    size_t result = 0;

    /** We open the "/proc/cpuinfo" file. */
    FILE* file = fopen ("/proc/cpuinfo", "r");
    if (file)
    {
        char buffer[256];
        while (fgets(buffer, sizeof(buffer), file))  {  if (strstr(buffer, "processor") != NULL)  { result ++;  }  }
        fclose (file);
    }

    if (result==0)  { result = 1; }

    return result;
}

/********************************************************************************/
string SystemInfoLinux::getHostName () const
{
    string result;

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname (hostname, sizeof(hostname)-1);
    result.assign (hostname, strlen(hostname));

    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryPhysicalTotal () const
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.totalram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryPhysicalUsed () const
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.totalram - memInfo.freeram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryBuffers () const
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.bufferram;
    result *= memInfo.mem_unit;
    return result;
}

#endif

/*********************************************************************
            #     #     #      #####   #######   #####
            ##   ##    # #    #     #  #     #  #     #
            # # # #   #   #   #        #     #  #
            #  #  #  #     #  #        #     #   #####
            #     #  #######  #        #     #        #
            #     #  #     #  #     #  #     #  #     #
            #     #  #     #   #####   #######   #####
*********************************************************************/

#ifdef __APPLE__

#include <sys/types.h>
#include <sys/sysctl.h>

#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>

/********************************************************************************/
size_t SystemInfoMacos::getNbCores () const
{
    int numCPU = 0;

    int mib[4];
    size_t len = sizeof(numCPU);

    /* set the mib for hw.ncpu */
    mib[0] = CTL_HW;
    mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;

    /* get the number of CPUs from the system */
    sysctl (mib, 2, &numCPU, &len, NULL, 0);

    if (numCPU < 1)
    {
        mib[1] = HW_NCPU;
        sysctl (mib, 2, &numCPU, &len, NULL, 0 );
    }

    if (numCPU < 1)  {  numCPU = 1;  }

    return numCPU;
}

/********************************************************************************/
string SystemInfoMacos::getHostName () const
{
    string result;

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname (hostname, sizeof(hostname)-1);
    result.assign (hostname, strlen(hostname));

    return result;
}

/********************************************************************************/
u_int64_t SystemInfoMacos::getMemoryPhysicalTotal () const
{
	int mib[2] = { CTL_HW, HW_MEMSIZE };
	size_t namelen = sizeof(mib) / sizeof(mib[0]);
	u_int64_t size;
	size_t len = sizeof(size);

	if (sysctl(mib, namelen, &size, &len, NULL, 0) < 0)
	{
		throw Exception ("Unable to get physical memory");
	}
	else
	{
		return size;
	}
}

/********************************************************************************/
u_int64_t SystemInfoMacos::getMemoryPhysicalUsed () const
{
	// see http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
	vm_size_t page_size;
	mach_port_t mach_port;
	mach_msg_type_number_t count;
	vm_statistics_data_t vm_stats;

	mach_port = mach_host_self();
	count = sizeof(vm_stats) / sizeof(natural_t);
	if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
	    KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO,
	                                    (host_info_t)&vm_stats, &count))
	{
		int64_t myFreeMemory = (int64_t)vm_stats.free_count * (int64_t)page_size;

		int64_t used_memory = ((int64_t)vm_stats.active_count +
	                   (int64_t)vm_stats.inactive_count +
	                   (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
		return myFreeMemory + used_memory;
	}
	else
	{
		throw Exception ("Unable to get free memory");
	}
}

#endif

/*********************************************************************
        #     #  ###  #     #  ######   #######  #     #   #####
        #  #  #   #   ##    #  #     #  #     #  #  #  #  #     #
        #  #  #   #   # #   #  #     #  #     #  #  #  #  #
        #  #  #   #   #  #  #  #     #  #     #  #  #  #   #####
        #  #  #   #   #   # #  #     #  #     #  #  #  #        #
        #  #  #   #   #    ##  #     #  #     #  #  #  #  #     #
         ## ##   ###  #     #  ######   #######   ## ##    #####
*********************************************************************/

#ifdef __WINDOWS__

#include <windows.h>

/********************************************************************************/
size_t SystemInfoWindows::getNbCores () const
{
    size_t result = 0;

    SYSTEM_INFO sysinfo;
    GetSystemInfo (&sysinfo);
    result = sysinfo.dwNumberOfProcessors;

    if (result==0)  { result = 1; }

    return result;
}

/********************************************************************************/
string SystemInfoWindows::getHostName () const
{
    string result;

    TCHAR  infoBuf[1024];
    DWORD  bufCharCount = sizeof(infoBuf)/sizeof(infoBuf[0]);

    if (GetComputerName( infoBuf, &bufCharCount ) )
    {
        result.assign (infoBuf, strlen(infoBuf));
    }

    return result;
}

#endif

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
