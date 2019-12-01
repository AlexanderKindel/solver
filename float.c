#include "declarations.h"

struct Float*float_copy(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_copy(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
    return out;
}

struct Float*float_reduce(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Integer*significand, struct Integer*exponent)
{
    if (!significand->value_count)
    {
        return &float_zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = ALLOCATE(output_stack, struct Float);
    while (!(significand->value[0] & 1))
    {
        significand = integer_halve(local_stack, significand);
        exponent = integer_add(local_stack, exponent, &one);
    }
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

struct Float*float_get_magnitude(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_get_magnitude(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
    return out;
}

struct Float*float_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b)
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
                integer_exponentiate(local_stack, output_stack, INT(2, 1),
                    exponent_difference)));
        out->exponent = integer_copy(output_stack, b->exponent);
    }
    else if (exponent_difference->sign < 0)
    {
        out = float_add(output_stack, local_stack, b, a);
    }
    else
    {
        out = float_reduce(output_stack, local_stack,
            integer_add(local_stack, a->significand, b->significand), a->exponent);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_negate(struct Stack*output_stack, struct Float*a)
{
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out->significand = integer_negate(output_stack, a->significand);
    out->exponent = integer_copy(output_stack, a->exponent);
    return out;
}

struct Float*float_subtract(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*minuend, struct Float*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = ALLOCATE(output_stack, struct Float);
    out = float_add(output_stack, local_stack, minuend, float_negate(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Float*float_multiply(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*out = float_reduce(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->significand, b->significand),
        integer_add(local_stack, a->exponent, b->exponent));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

int8_t float_get_sign(struct Float*a)
{
    return a->significand->sign;
}

int8_t float_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Float*a, struct Float*b)
{
    void*local_stack_a_savepoint = local_stack_a->cursor;
    int8_t out = float_subtract(local_stack_a, local_stack_b, a, b)->significand->sign;
	local_stack_a->cursor = local_stack_a_savepoint;
    return out;
}

struct Float*float_get_min(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b)
{
    if (float_compare(output_stack, local_stack, a, b) < 0)
    {
        return float_copy(output_stack, a);
    }
    else
    {
        return float_copy(output_stack, b);
    }
}

struct Float*float_get_max(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float*a, struct Float*b)
{
    if (float_compare(output_stack, local_stack, a, b) > 0)
    {
        return float_copy(output_stack, a);
    }
    else
    {
        return float_copy(output_stack, b);
    }
}

void float_estimate_root(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Float**out_get_min, struct Float**out_get_max, struct Float*a, struct Rational*interval_size,
    struct Integer*index)
{
    ASSERT(a->significand->sign >= 0, "float_estimate_root was called on a negative a value.");
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*rational_radicand = float_to_rational(local_stack, output_stack, a);
    if (rational_compare(output_stack, local_stack, rational_radicand, &rational_one) < 0)
    {
        *out_get_max = &float_one;
    }
    else
    {
        *out_get_max = float_copy(local_stack, a);
    }
    struct Integer*index_get_minus_one = integer_add(local_stack, index, INT(1, -1));
    while (true)
    {
        struct Rational*delta = rational_integer_divide(local_stack, output_stack,
            rational_subtract(local_stack, output_stack,
                rational_divide(local_stack, output_stack, rational_radicand,
                    float_to_rational(local_stack, output_stack,
                        float_exponentiate(local_stack, output_stack, *out_get_max, index_get_minus_one))),
                float_to_rational(local_stack, output_stack, *out_get_max)), index);
        struct FloatInterval*delta_float_estimate =
            rational_get_float_estimate(local_stack, output_stack, delta,
                rational_integer_divide(local_stack, output_stack, delta, INT(2, -1)));
        *out_get_max = float_add(local_stack, output_stack, *out_get_max, delta_float_estimate->max);
        if (rational_compare(output_stack, local_stack,
            rational_subtract(local_stack, output_stack,
                float_to_rational(local_stack, output_stack, delta_float_estimate->max),
                rational_double(local_stack, output_stack, delta)), interval_size) <= 0)
        {
            *out_get_min = float_add(output_stack, local_stack, *out_get_max, delta_float_estimate->max);
            *out_get_max = float_copy(output_stack, *out_get_max);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
    }
}

struct Rational*float_to_rational(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Float*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    if (a->exponent->sign < 0)
    {
        out->numerator = integer_copy(output_stack, a->significand);
        out->denominator = integer_exponentiate(output_stack, local_stack, INT(2, 1),
            integer_get_magnitude(local_stack, a->exponent));
    }
    else
    {
        out->numerator = integer_multiply(output_stack, local_stack,
            integer_exponentiate(local_stack, output_stack, INT(2, 1), a->exponent),
            a->significand);
        out->denominator = &one;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}