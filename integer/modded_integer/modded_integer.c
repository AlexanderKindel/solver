#include "declarations.h"

struct Integer*integer_get_residue(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out =
        integer_euclidean_divide(local_stack, output_stack, a, characteristic).remainder;
    if (out->sign < 0)
    {
        out = integer_add(output_stack, out, characteristic);
    }
    else
    {
        out = integer_copy(output_stack, out);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*modded_integer_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = integer_euclidean_divide(output_stack, local_stack,
        integer_add(local_stack, a, b), characteristic).remainder;
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*modded_integer_negate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out =
        integer_get_residue(output_stack, local_stack, integer_negate(local_stack, a), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*modded_integer_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*minuend, struct Integer*subtrahend,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = integer_get_residue(output_stack, local_stack,
        integer_subtract(local_stack, output_stack, minuend, subtrahend), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*modded_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = integer_euclidean_divide(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a, b), characteristic).remainder;
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Integer*modded_integer_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo info =
        integer_get_extended_gcd(local_stack, output_stack, a, characteristic);
    struct Integer*out;
    if (((struct Integer*)info.a_coefficient)->sign < 0)
    {
        out = integer_add(output_stack, info.a_coefficient, characteristic);
    }
    else
    {
        out = integer_copy(output_stack, info.a_coefficient);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

#include "integer/modded_integer/modded_polynomial/modded_polynomial.c"