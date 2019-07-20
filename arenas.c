#include "declarations.h"

size_t next_page_boundary(size_t address)
{
    size_t distance_from_page_start = address % page_size;
    if (distance_from_page_start)
    {
        return address + page_size - distance_from_page_start;
    }
    return address;
}

size_t next_aligned_address(size_t address, size_t alignment)
{
    return ((address + alignment - 1) / alignment) * alignment;
}

void stack_initialize(struct Stack*out, size_t start, size_t end)
{
    out->start = start;
    out->end = end;
    out->cursor = (void*)start;
    out->allocation_cursor = out->start;
}

void stack_free(struct Stack*out)
{
    while (out->allocation_cursor > out->start)
    {
        out->allocation_cursor -= page_size;
        VirtualFree((void*)out->allocation_cursor, 0, MEM_RELEASE);
    }
    out->cursor = (void*)out->start;
}

void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment)
{
    size_t slot = next_aligned_address((size_t)output_stack->cursor, alignment);
    output_stack->cursor = (void*)(slot + slot_size);
    if ((size_t)output_stack->cursor > output_stack->end)
    {
        crash("Stack ran out of virtual address space.");
    }
    while ((size_t)output_stack->cursor > output_stack->allocation_cursor)
    {
        if (!VirtualAlloc((void*)output_stack->allocation_cursor, page_size,
            MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE))
        {
            crash("Ran out of physical memory.");
        }
        output_stack->allocation_cursor = output_stack->allocation_cursor + page_size;
    }
    return (void*)slot;
}

void pool_set_initialize(struct PoolSet*out, size_t start, size_t end)
{
    out->pool_page =
        VirtualAlloc((void*)start, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    stack_initialize(&out->slot_page_stack, start + page_size, end);
}

void*pool_slot_page_allocate(struct PoolSet*pool_set)
{
    pool_set->slot_page_stack.allocation_cursor =
        pool_set->slot_page_stack.allocation_cursor + page_size;
    if (pool_set->slot_page_stack.allocation_cursor > pool_set->slot_page_stack.end)
    {
        crash("Pool stack ran out of virtual address space.");
    }
    void*new_page_address = pool_set->slot_page_stack.cursor;
    pool_set->slot_page_stack.cursor = (void*)pool_set->slot_page_stack.allocation_cursor;
    return VirtualAlloc(new_page_address, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
}

struct Pool*get_pool(struct PoolSet*pool_set, size_t slot_size)
{
    return pool_set->pool_page + (slot_size + sizeof(size_t) - 1) / sizeof(size_t) - 2;
}

//Produces pointer-aligned memory blocks, which would be incorrect for allocating types whose
//alignment is larger than a pointer, but the program currently uses no such types.
void*pool_slot_allocate(struct PoolSet*pool_set, size_t slot_size)
{
    struct Pool*pool = get_pool(pool_set, slot_size);
    if (!pool->first_page)
    {
        pool->first_page = pool_slot_page_allocate(pool_set);
        pool->cursor = (void*)((size_t)pool->first_page + sizeof(void*));
    }
    void*out = pool->free_list;
    if (out)
    {
        pool->free_list = *(void**)out;
    }
    else
    {
        out = pool->cursor;
        pool->cursor = (void*)((size_t)pool->cursor + slot_size);
        size_t boundary = next_page_boundary((size_t)out);
        if ((size_t)pool->cursor > boundary)
        {
            void*new_page = pool_slot_page_allocate(pool_set);
            *(void**)(boundary - page_size) = new_page;
            out = (void*)((size_t)new_page + sizeof(void*));
        }
    }
    return out;
}

void pool_slot_free(struct PoolSet*pool_set, void*a, size_t slot_size)
{
    ASSERT(pool_memory_start <= (size_t)a && (size_t)a <= pool_memory_end,
        "pool_slot_free argument wasn't in a pool.\n");
    struct Pool*pool = get_pool(pool_set, slot_size);
    *(void**)a = pool->free_list;
    pool->free_list = a;
}