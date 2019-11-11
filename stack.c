#include "declarations.h"

__declspec(noreturn) void crash(char*message)
{
    puts(message);
    abort();
}

void stack_initialize(struct Stack*out, size_t start, size_t size)
{
    out->start = start;
    out->end = start + size;
    out->cursor = (void*)start;
    out->cursor_max = out->start;
}

void stack_free(struct Stack*out)
{
    DECOMMIT_STACK(out);
    out->cursor = (void*)out->start;
    out->cursor_max = out->start;
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
        puts("Insufficient memory allocated for calculation.");
        longjmp(memory_error_buffer, 0);
    }
    COMMIT_MEMORY(output_stack);
}

void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment)
{
    void*slot = array_start(output_stack, alignment);
    extend_array(output_stack, slot_size);
    return slot;
}