#include "declarations.h"

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