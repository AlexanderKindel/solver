#ifndef OS_H
#define OS_H

size_t page_size;

#if defined(_WIN64)

#define restrict __restrict

#include <windows.h>
#include <intrin.h>

//Ensures that VirtualAlloc never rounds addresses that are multiples of page_size.
#define SET_PAGE_SIZE()\
{\
    SYSTEM_INFO system_info;\
    GetSystemInfo(&system_info);\
    page_size = max(system_info.dwAllocationGranularity, system_info.dwPageSize);\
}

#define RESERVE_MEMORY(size) VirtualAlloc(0, size, MEM_RESERVE, PAGE_READWRITE)

#define COMMIT_PAGE(address) VirtualAlloc(address, page_size, MEM_COMMIT, PAGE_READWRITE)

#define DECOMMIT_STACK(stack)\
VirtualFree(stack->start, (uintptr_t)stack->cursor_max - (uintptr_t)stack->start, MEM_DECOMMIT)

uint32_t get_leading_digit_place(uint32_t a)
{
    uint32_t out;
    _BitScanReverse64(&out, a);
    return out;
}

#elif __has_include (<unistd.h>)

#include <sys/mman.h>
#include <unistd.h>

#define SET_PAGE_SIZE() page_size = sysconf(_SC_PAGESIZE)

#define RESERVE_MEMORY(size) mmap(0, size, PROT_NONE, MAP_ANON | MAP_PRIVATE, -1, 0)

#define COMMIT_PAGE(address) mprotect(address, page_size, PROT_READ | PROT_WRITE)

#define DECOMMIT_STACK(stack)\
mprotect(stack->start, (uintptr_t)stack->cursor_max - (uintptr_t)stack->start, PROT_NONE)

uint32_t get_leading_digit_place(uint32_t a)
{
    uint32_t divisor_leading_digit = 0x80000000;
    uint32_t i = 31;
    while (true)
    {
        if ((a & divisor_leading_digit) != 0)
        {
            return i;
        }
        divisor_leading_digit = divisor_leading_digit >> 1;
        --i;
    }
}

#else

#error unsupported platform

#endif

#endif
