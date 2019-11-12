#ifndef OS_H
#define OS_H

size_t page_size;

#ifdef _WIN64

#include <windows.h>

//Ensures that VirtualAlloc never rounds addresses that are multiples of page_size.
#define SET_PAGE_SIZE()\
{\
    SYSTEM_INFO system_info;\
    GetSystemInfo(&system_info);\
    page_size = max(system_info.dwAllocationGranularity, system_info.dwPageSize);\
}

#define RESERVE_MEMORY(size) VirtualAlloc(0, size, MEM_RESERVE, PAGE_READWRITE)

#define COMMIT_MEMORY(stack)\
while ((size_t)stack->cursor > stack->cursor_max)\
{\
    if (!VirtualAlloc((void*)stack->cursor_max, page_size, MEM_COMMIT, PAGE_READWRITE))\
    {\
        puts("Insufficient physical memory for calculation.");\
        longjmp(memory_error_buffer, 0);\
    }\
    stack->cursor_max = stack->cursor_max + page_size;\
}

#define DECOMMIT_STACK(stack)\
VirtualFree((void*)stack->start, stack->cursor_max - stack->start, MEM_DECOMMIT)

#endif

#endif