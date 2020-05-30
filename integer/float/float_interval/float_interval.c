#include "declarations.h"

struct FloatInterval float_interval_copy(struct Stack*output_stack, struct FloatInterval*a)
{
    return (struct FloatInterval) { float_copy(output_stack, &a->min),
        float_copy(output_stack, &a->max) };
}

struct Float float_interval_get_max_magnitude(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float min_magnitude = float_get_magnitude(local_stack, &a->min);
    struct Float max_magnitude = float_get_magnitude(local_stack, &a->max);
    struct Float out = float_get_max(output_stack, local_stack, &min_magnitude, &max_magnitude);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval float_interval_add(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b)
{
    return (struct FloatInterval) { float_add(output_stack, local_stack, &a->min, &b->min),
        float_add(output_stack, local_stack, &a->max, &b->max) };
}

struct FloatInterval float_interval_negate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a)
{
    struct FloatInterval out;
    if (float_get_sign(&a->min) != -float_get_sign(&a->max))
    {
        out.min = float_negate(output_stack, &a->max);
        out.max = float_negate(output_stack, &a->min);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Float max_negation = float_negate(local_stack, &a->max);
        out.min = float_get_min(output_stack, local_stack, &max_negation, &a->min);
        struct Float min_negation = float_negate(local_stack, &a->min);
        out.max = float_get_max(output_stack, local_stack, &min_negation, &a->max);
        local_stack->cursor = local_stack_savepoint;
    }
    return out;
}

struct FloatInterval float_interval_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct FloatInterval*minuend,
    struct FloatInterval*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval negative_subtrahend =
        float_interval_negate(local_stack, output_stack, subtrahend);
    struct FloatInterval out =
        float_interval_add(output_stack, local_stack, minuend, &negative_subtrahend);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval float_interval_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b)
{
    struct FloatInterval out;
    if (float_get_sign(&a->min) >= 0)
    {
        if (float_get_sign(&b->min) >= 0)
        {
            out.min = float_multiply(output_stack, local_stack, &a->min, &b->min);
            out.max = float_multiply(output_stack, local_stack, &a->max, &b->max);
        }
        else
        {
            out.min = float_multiply(output_stack, local_stack, &a->max, &b->min);
            if (float_get_sign(&b->max) <= 0)
            {
                out.max = float_multiply(output_stack, local_stack, &a->min, &b->max);
            }
            else
            {
                out.max = float_multiply(output_stack, local_stack, &a->max, &b->max);
            }
        }
    }
    else if (float_get_sign(&a->max) <= 0)
    {
        if (float_get_sign(&b->min) >= 0)
        {
            out.min = float_multiply(output_stack, local_stack, &a->min, &b->max);
            out.max = float_multiply(output_stack, local_stack, &a->max, &b->min);
        }
        else
        {
            if (float_get_sign(&b->max) <= 0)
            {
                out.min = float_multiply(output_stack, local_stack, &a->max, &b->max);
            }
            else
            {
                out.min = float_multiply(output_stack, local_stack, &a->min, &b->max);
            }
            out.max = float_multiply(output_stack, local_stack, &a->min, &b->min);
        }
    }
    else
    {
        if (float_get_sign(&b->min) >= 0)
        {
            out.min = float_multiply(output_stack, local_stack, &a->min, &b->max);
            out.max = float_multiply(output_stack, local_stack, &a->max, &b->max);
        }
        else if (float_get_sign(&b->max) <= 0)
        {
            out.min = float_multiply(output_stack, local_stack, &a->max, &b->min);
            out.max = float_multiply(output_stack, local_stack, &a->min, &b->min);
        }
        else
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Float product_a = float_multiply(local_stack, output_stack, &a->min, &b->max);
            struct Float product_b = float_multiply(local_stack, output_stack, &a->max, &b->min);
            out.min = float_get_min(output_stack, local_stack, &product_a, &product_b);
            product_a = float_multiply(local_stack, output_stack, &a->min, &b->min);
            product_b = float_multiply(local_stack, output_stack, &a->max, &b->max);
            out.max = float_get_max(output_stack, local_stack, &product_a, &product_b);
            local_stack->cursor = local_stack_savepoint;
        }
    }
    return out;
}

bool float_intervals_are_disjoint(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct FloatInterval*a, struct FloatInterval*b)
{
    return float_compare(local_stack_a, local_stack_b, &a->min, &b->max) > 0 ||
        float_compare(local_stack_a, local_stack_b, &b->min, &a->max) > 0;
}

bool float_interval_a_contains_b(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct FloatInterval*a, struct FloatInterval*b)
{
    return float_compare(local_stack_a, local_stack_b, &a->min, &b->min) <= 0 &&
        float_compare(local_stack_a, local_stack_b, &a->max, &b->max) >= 0;
}

struct RationalInterval float_interval_to_rational_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct FloatInterval*a)
{
    return (struct RationalInterval) { float_to_rational(output_stack, local_stack, &a->min),
        float_to_rational(output_stack, local_stack, &a->max) };
}

#include "integer/float/float_interval/region/region.c"