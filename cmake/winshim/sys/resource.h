/**
 * part of LukasNalbach/Move-r
 *
 * Windows shim for <sys/resource.h>. Only reachable on Windows builds (this directory is on
 * the include path solely inside `if(WIN32)`). MinGW's libstdc++/POSIX build of sdsl includes
 * this header for getrusage()/ru_maxrss; Windows has no getrusage, so this stub reports 0
 * (peak-memory statistics are unavailable, which matches the malloc_count stub behaviour).
 */
#pragma once

#if defined(_WIN32)

#include <cstring>

#define RUSAGE_SELF     0
#define RUSAGE_CHILDREN (-1)

struct rusage {
    struct { long tv_sec; long tv_usec; } ru_utime;
    struct { long tv_sec; long tv_usec; } ru_stime;
    long ru_maxrss;
    long ru_ixrss;
    long ru_idrss;
    long ru_isrss;
    long ru_minflt;
    long ru_majflt;
    long ru_nswap;
    long ru_inblock;
    long ru_oublock;
    long ru_msgsnd;
    long ru_msgrcv;
    long ru_nsignals;
    long ru_nvcsw;
    long ru_nivcsw;
};

static inline int getrusage(int, struct rusage* r)
{
    if (r) std::memset(r, 0, sizeof(*r));
    return 0;
}

#endif // _WIN32
