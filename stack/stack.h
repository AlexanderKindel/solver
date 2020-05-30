/*Let a function's result refer to its return value together with the values to which it assigns its
out parameters, if it has any. Let the reference graph of a value of a primitive data type refer to
the value, the reference graph of a pointer refer to the reference graph of the value it points to,
and the reference graph of an aggregate data type, to the union of the reference graphs of its
elements.

If the reference graph of a function's result contains dynamically allocated elements, then the
function takes an output_stack parameter, and allocates all such elements on that Stack. This means
that functions never assign dynamically allocated values passed in by their callers directly into
the reference graphs of their results, instead deep-copying such values to output_stack first.

Functions that need to do dynamic allocations in service of calculating their results but that are
not themselves part of the result occasionally do such allocations on output_stack for efficiency
reasons, but most of the time, they take a local_stack parameter and do the allocations on that. Any
allocations a function does on a local_stack are popped before the function returns. The few
functions that do a significant amount of allocation on output_stack that isn't part of their
results are identified with comments.*/

void stack_initialize(struct Stack*out, void*start, size_t size);
void stack_reset(struct Stack*stack);
void*align_cursor(struct Stack*output_stack, size_t alignment);
void unaligned_stack_slot_allocate(struct Stack*output_stack, size_t element_size);
void*stack_slot_allocate(struct Stack*output_stack, size_t slot_size, size_t alignment);