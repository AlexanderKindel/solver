#include "declarations.h"

__declspec(noreturn) void crash(char*message)
{
    puts(message);
    abort();
}

void stack_initialize(struct Stack*out, void*start, size_t size)
{
    out->start = start;
    out->end = (uint8_t*)start + size;
    out->cursor = start;
    out->cursor_max = out->start;
}

void stack_reset(struct Stack*stack)
{
    DECOMMIT_STACK(stack);
    stack->cursor = stack->start;
    stack->cursor_max = stack->start;
}

void*array_start(struct Stack*output_stack, size_t alignment)
{
    output_stack->cursor =
        (void*)(((uintptr_t)output_stack->cursor + alignment - 1) & -(uintptr_t)alignment);
    return output_stack->cursor;
}

void array_extend(struct Stack*output_stack, size_t element_size)
{
    output_stack->cursor = (uint8_t*)output_stack->cursor + element_size;
    if (output_stack->cursor > output_stack->end)
    {
        puts("Insufficient memory allocated for calculation.");
        longjmp(memory_error_buffer, 0);
    }
    while (output_stack->cursor > output_stack->cursor_max)
    {
        COMMIT_PAGE(output_stack->cursor_max);
        output_stack->cursor_max = (uint8_t*)output_stack->cursor_max + page_size;
    }
}

void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment)
{
    void*slot = array_start(output_stack, alignment);
    array_extend(output_stack, slot_size);
    return slot;
}