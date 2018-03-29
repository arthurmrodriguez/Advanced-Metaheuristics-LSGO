#ifndef FREEMEM_H_
#define FREEMEM_H_

#include "SingleLogStat.h"

#if __APPLE__

#import <sys/sysctl.h>
#import <mach/host_info.h>
#import <mach/mach_host.h>
#import <mach/task_info.h>
#import <mach/task.h>

#elif __linux

#include <sys/sysinfo.h>

#endif



class FreememLogStat : public SingleLogStat {

public:
	FreememLogStat(const Algorithm& alg) : SingleLogStat("freemem",alg) {}
	~FreememLogStat(){}
	
  long double computeValue(){

#if __APPLE__
    // From http://stackoverflow.com/questions/6094444/how-can-i-programmatically-check-free-system-memory-on-mac-like-the-activity-mon
    vm_statistics_data_t   vmstat;
    mach_msg_type_number_t count = HOST_VM_INFO_COUNT;

    if (host_statistics (mach_host_self (), HOST_VM_INFO, (host_info_t) &vmstat, &count) != KERN_SUCCESS) {
        fprintf (stderr, "Failed to get VM statistics.");
    }

    long double total = vmstat.wire_count + vmstat.active_count + vmstat.inactive_count + vmstat.free_count;
    long double free  = vmstat.free_count / total;

    return free;

#elif __linux
    struct sysinfo info;



    if (sysinfo(&info) != 0) throw runtime_error("sysinfo: error reading system statistics");

    return info.freeram;
#endif

  }
};

#endif /*FREEMEM_H_*/
