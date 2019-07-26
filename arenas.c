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

void stack_initialize(struct Stack*out, size_t start, size_t end)
{
    out->start = start;
    out->end = end;
    out->cursor = (void*)start;
    out->cursor_max = out->start;
}

void stack_free(struct Stack*out)
{
    while (out->cursor_max > out->start)
    {
        out->cursor_max -= page_size;
        VirtualFree((void*)out->cursor_max, 0, MEM_RELEASE);
    }
    out->cursor = (void*)out->start;
}

void*array_start(struct Stack*output_stack, size_t alignment)
{
    output_stack->cursor =
        (void*)((((size_t)output_stack->cursor + alignment - 1) / alignment) * alignment);
    return output_stack->cursor;
}

void extend_array(struct Stack*output_stack, size_t element_size)
{
    output_stack->cursor = (void*)((size_t)output_stack->cursor + element_size);
    if ((size_t)output_stack->cursor > output_stack->end)
    {
        crash("Stack ran out of virtual address space.");
    }
    while ((size_t)output_stack->cursor > output_stack->cursor_max)
    {
        if (!VirtualAlloc((void*)output_stack->cursor_max, page_size, MEM_RESERVE | MEM_COMMIT,
            PAGE_READWRITE))
        {
            crash("Ran out of physical memory.");
        }
        output_stack->cursor_max = output_stack->cursor_max + page_size;
    }
}

void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment)
{
    void*slot = array_start(output_stack, alignment);
    extend_array(output_stack, slot_size);
    return slot;
}

void pool_set_initialize(struct PoolSet*out, size_t start, size_t end)
{
    out->pool_page =
        VirtualAlloc((void*)start, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    stack_initialize(&out->slot_page_stack, start + page_size, end);
}

struct PoolSlotPage*pool_slot_page_allocate(struct PoolSet*pool_set)
{
    pool_set->slot_page_stack.cursor_max = pool_set->slot_page_stack.cursor_max + page_size;
    if (pool_set->slot_page_stack.cursor_max > pool_set->slot_page_stack.end)
    {
        crash("Pool stack ran out of virtual address space.");
    }
    void*new_page_address = pool_set->slot_page_stack.cursor;
    pool_set->slot_page_stack.cursor = (void*)pool_set->slot_page_stack.cursor_max;
    return VirtualAlloc(new_page_address, page_size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
}

struct Pool*get_pool(struct PoolSet*pool_set, size_t value_size)
{
    return pool_set->pool_page + (value_size + sizeof(size_t) - 1) / sizeof(size_t) - 2;
}

//Produces pointer-aligned values, which would be incorrect for allocating types whose alignments
//are larger than a pointer, but the program currently uses no such types.
void*pool_value_allocate(struct PoolSet*pool_set, size_t size)
{
    struct Pool*pool = get_pool(pool_set, size);
    if (!pool->first_page)
    {
        pool->first_page = pool_slot_page_allocate(pool_set);
        pool->cursor = &pool->first_page->first_slot;
    }
    struct PoolSlot*out = pool->free_list;
    if (out)
    {
        pool->free_list = out->next;
    }
    else
    {
        out = pool->cursor;
        pool->cursor =
            (struct PoolSlot*)((size_t)pool->cursor + offsetof(struct PoolSlot, value) + size);
        size_t boundary = next_page_boundary((size_t)out);
        if ((size_t)pool->cursor > boundary)
        {
            struct PoolSlotPage*new_page = pool_slot_page_allocate(pool_set);
            ((struct PoolSlotPage*)(boundary - page_size))->next = new_page;
            out = &new_page->first_slot;
        }
    }
    out->reference_count = 1;
    return out->value;
}

struct PoolSlot*pool_slot_from_value(void*a)
{
    return (struct PoolSlot*)((size_t)a - offsetof(struct PoolSlot, value));
}

void increment_reference_count(void*a)
{
    ASSERT(pool_memory_start <= (size_t)a && (size_t)a <= pool_memory_end,
        "increment_reference_count argument wasn't in a pool.\n");
    struct PoolSlot*slot = pool_slot_from_value(a);
    ++slot->reference_count;
}

void increment_reference_count_if_non_null(void*a)
{
    if (a)
    {
        increment_reference_count(a);
    }
}

void pool_slot_free(struct PoolSet*pool_set, struct PoolSlot*a, size_t size)
{
    ASSERT(pool_memory_start <= (size_t)a && (size_t)a <= pool_memory_end,
        "pool_slot_free argument wasn't in a pool.\n");
    ASSERT(!a->reference_count, "pool_slot_free called on slot with nonzero reference count.");
    struct Pool*pool = get_pool(pool_set, size);
    a->next = pool->free_list;
    pool->free_list = a;
}