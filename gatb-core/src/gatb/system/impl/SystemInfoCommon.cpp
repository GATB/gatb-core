/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Copyright (c) 2013                                                      *
 *                                                                           *
 *   GATB is free software; you can redistribute it and/or modify it under   *
 *   the CECILL version 2 License, that is compatible with the GNU General   *
 *   Public License                                                          *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   CECILL version 2 License for more details.                              *
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
size_t SystemInfoLinux::getNbCores ()
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
string SystemInfoLinux::getHostName ()
{
    string result;

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname (hostname, sizeof(hostname)-1);
    result.assign (hostname, strlen(hostname));

    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryPhysicalTotal ()
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.totalram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryPhysicalUsed ()
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.totalram - memInfo.freeram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryBuffers ()
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

/********************************************************************************/
size_t SystemInfoMacos::getNbCores ()
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
string SystemInfoMacos::getHostName ()
{
    string result;

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname (hostname, sizeof(hostname)-1);
    result.assign (hostname, strlen(hostname));

    return result;
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
size_t SystemInfoWindows::getNbCores ()
{
    size_t result = 0;

    SYSTEM_INFO sysinfo;
    GetSystemInfo (&sysinfo);
    result = sysinfo.dwNumberOfProcessors;

    if (result==0)  { result = 1; }

    return result;
}

/********************************************************************************/
string SystemInfoWindows::getHostName ()
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
