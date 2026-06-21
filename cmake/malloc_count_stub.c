/******************************************************************************
 * malloc_count_stub.c
 *
 * No-op replacement for external/malloc_count/malloc_count.c, compiled in place
 * of it when the CMake option MOVE_R_USE_MALLOC_COUNT is OFF. It provides the
 * same public symbols declared in malloc_count.h, but does not override
 * malloc()/free(), so the allocation-tracking overhead is avoided entirely.
 *
 * All measurement functions return 0; memory statistics reported by move-r are
 * therefore reported as 0 B when malloc_count is disabled.
 *****************************************************************************/

#include <malloc_count.h>

extern size_t malloc_count_current(void)
{
    return 0;
}

extern size_t malloc_count_peak(void)
{
    return 0;
}

extern void malloc_count_reset_peak(void)
{
}

extern size_t malloc_count_num_allocs(void)
{
    return 0;
}

void malloc_count_set_callback(malloc_count_callback_type cb, void* cookie)
{
    (void)cb;
    (void)cookie;
}

extern void malloc_count_print_status(void)
{
}
