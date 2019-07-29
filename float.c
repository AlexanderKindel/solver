#include "declarations.h"

struct Float*float_copy(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_copy(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
    return out;
}

struct Float*float_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*significand, struct Integer*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!(significand->value[0] & 1))
    {
        significand = integer_half(local_stack, significand);
        exponent = integer_add(local_stack, exponent, &one);
    }
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_copy(output_stack, significand);
    out->exponent = integer_copy(output_stack, exponent);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

bool float_equals(struct Float*a, struct Float*b)
{
    return integer_equals(a->significand, b->significand) &&
        integer_equals(a->exponent, b->exponent);
}

struct Float*float_magnitude(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_magnitude(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
    return out;
}

struct Float*float_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*exponent_difference =
        integer_subtract(local_stack, output_stack, a->exponent, b->exponent);
    struct Float*out;
    if (exponent_difference->sign > 0)
    {
        out = ALLOCATE(output_stack, struct Float);
        out->significand = integer_add(output_stack, b->significand,
            integer_multiply(local_stack, output_stack, a->significand,
                integer_exponentiate(local_stack, output_stack, &INT(2, +),
                    exponent_difference)));
        out->exponent = integer_copy(output_stack, b->exponent);
    }
    else if (exponent_difference->sign < 0)
    {
        out = float_add(output_stack, local_stack, b, a);
    }
    else
    {
        out = float_reduced(output_stack, local_stack,
            integer_add(local_stack, a->significand, b->significand), a->exponent);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_generic_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b, void*unused)
{
    return float_add(output_stack, local_stack, a, b);
}

struct Float*float_negative(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_negative(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
    return out;
}

struct Float*float_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    struct Float*a, void*unused)
{
    return float_negative(output_stack, a);
}

struct Float*float_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*minuend, struct Float*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out = float_add(output_stack, local_stack, minuend, float_negative(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = float_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->significand, b->significand),
        integer_add(local_stack, a->exponent, b->exponent));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a, struct Float*b, void*unused)
{
    return float_multiply(output_stack, local_stack, a, b);
}

int8_t float_compare(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t out = float_subtract(stack_a, stack_b, a, b)->significand->sign;
    stack_a->cursor = stack_a_savepoint;
    return out;
}

struct Float*float_max(struct Stack*stack_a, struct Stack*stack_b, struct Float*a, struct Float*b)
{
    if (float_compare(stack_a, stack_b, a, b) > 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

struct Float*float_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*base, struct Integer*exponent)
{
    return generic_exponentiate(&float_operations, output_stack, local_stack, base, exponent, 0);
}

void float_estimate_root(struct Stack*output_stack, struct Stack*local_stack, struct Float**out_min,
    struct Float**out_max, struct Float*a, struct Rational*interval_size, struct Integer*index)
{
    ASSERT(a->significand->sign >= 0, "float_estimate_root was called on a negative a value.\n");
    void*local_stack_savepoint = local_stack->cursor;
    if (a->exponent < 0)
    {
        (*out_max)->significand = &one;
        (*out_max)->exponent = &zero;
    }
    else
    {
        *out_max = float_copy(local_stack, a);
    }
    struct Rational*rational_radicand = float_to_rational(local_stack, output_stack, a);
    struct Integer*index_minus_one = integer_add(local_stack, index, &INT(1, -));
    while (true)
    {
        struct Rational*delta = rational_integer_divide(local_stack, output_stack,
            rational_subtract(local_stack, output_stack,
                rational_divide(local_stack, output_stack, rational_radicand,
                    float_to_rational(local_stack, output_stack,
                        float_exponentiate(local_stack, output_stack, *out_max, index_minus_one))),
                float_to_rational(local_stack, output_stack, *out_max)), index);
        struct FloatInterval delta_float_estimate;
        rational_float_estimate(local_stack, output_stack, &delta_float_estimate, delta,
            rational_integer_divide(local_stack, output_stack, delta, &INT(2, -)));
        *out_max = float_add(local_stack, output_stack, *out_max, delta_float_estimate.max);
        if (rational_compare(output_stack, local_stack,
            rational_subtract(local_stack, output_stack,
                float_to_rational(local_stack, output_stack, delta_float_estimate.max),
                rational_doubled(local_stack, output_stack, delta)), interval_size) <= 0)
        {
            *out_min = float_add(output_stack, local_stack, *out_max, delta_float_estimate.max);
            *out_max = float_copy(output_stack, *out_max);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
    }
}

struct Rational*float_to_rational(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    if (a->exponent->sign < 0)
    {
        out->numerator = integer_copy(output_stack, a->significand);
        out->denominator = integer_exponentiate(output_stack, local_stack, &INT(2, +),
            integer_magnitude(local_stack, a->exponent));
    }
    else
    {
        out->numerator = integer_multiply(output_stack, local_stack,
            integer_exponentiate(local_stack, output_stack, &INT(2, +), a->exponent),
            a->significand);
        out->denominator = &one;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

bool float_intervals_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval*a, struct FloatInterval*b)
{
    return float_compare(stack_a, stack_b, a->min, b->max) > 0 ||
        float_compare(stack_a, stack_b, b->min, a->max) > 0;
}

void float_interval_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct FloatInterval*a, struct FloatInterval*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->min = float_multiply(local_stack, output_stack, a->min, b->min);
    out->max = out->min;
    struct Float*bound_candidates[3] =
    { float_multiply(local_stack, output_stack, a->min, b->max),
        float_multiply(local_stack, output_stack, a->max, b->min),
        float_multiply(local_stack, output_stack, a->max, b->max) };
    for (uint8_t i = 0; i < 3; ++i)
    {
        if (float_compare(output_stack, local_stack, bound_candidates[i], out->min) <= 0)
        {
            out->min = bound_candidates[i];
        }
        else if (float_compare(output_stack, local_stack, out->max, bound_candidates[i]) < 0)
        {
            out->max = bound_candidates[i];
        }
    }
    out->min = float_copy(output_stack, out->min);
    out->max = float_copy(output_stack, out->max);
    local_stack->cursor = local_stack_savepoint;
}

void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct FloatInterval*a)
{
    out->min = float_to_rational(output_stack, local_stack, a->min);
    out->max = float_to_rational(output_stack, local_stack, a->max);
}