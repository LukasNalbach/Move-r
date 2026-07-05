/**
 * part of LukasNalbach/Move-r
 *
 * Windows shim for <sys/mman.h>. Only reachable on Windows builds (this directory is
 * added to the include path solely inside `if(WIN32)` in the top-level CMakeLists, so it
 * never shadows the real header on Linux/macOS). It maps the small mmap-based growable
 * allocator API used by the bundled sux library (sux/util/Vector.hpp, Expandable.hpp)
 * onto calloc/realloc/free. This drops the transparent-hugepage optimization but is
 * functionally correct; sux is the only consumer.
 */
#pragma once

#if defined(_WIN32)

#include <cstddef>
#include <cstdlib>

#define PROT_NONE 0
#define PROT_READ 1
#define PROT_WRITE 2
#define MAP_SHARED 0x01
#define MAP_PRIVATE 0x02
#define MAP_ANONYMOUS 0x20
#define MAP_ANON MAP_ANONYMOUS
#define MAP_NORESERVE 0
#define MAP_HUGETLB 0
#define MAP_FAILED ((void*) -1)
#define MADV_NORMAL 0
#define MADV_HUGEPAGE 0
#define MREMAP_MAYMOVE 1

// mmap of anonymous memory: return zero-initialized storage (as MAP_ANONYMOUS would).
static inline void* mmap(void*, size_t length, int, int, int, long)
{
    void* p = std::calloc(1, length);
    return p ? p : MAP_FAILED;
}

static inline int munmap(void* addr, size_t)
{
    std::free(addr);
    return 0;
}

static inline void* mremap(void* old_address, size_t, size_t new_size, int)
{
    return std::realloc(old_address, new_size);
}

static inline int madvise(void*, size_t, int) { return 0; }

static inline int mprotect(void*, size_t, int) { return 0; }

#endif // _WIN32
