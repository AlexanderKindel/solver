#include "declarations.h"

struct RationalInterval rational_interval_copy(struct Stack*output_stack, struct RationalInterval*a)
{
    return (struct RationalInterval) { rational_copy(output_stack, &a->min),
        rational_copy(output_stack, &a->max) };
}

struct Rational rational_interval_get_max_magnitude(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational min_magnitude = rational_get_magnitude(local_stack, &a->min);
    struct Rational max_magnitude = rational_get_magnitude(local_stack, &a->max);
    struct Rational out =
        rational_get_max(output_stack, local_stack, &min_magnitude, &max_magnitude);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval rational_interval_add(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*a, struct RationalInterval*b)
{
    return (struct RationalInterval) { rational_add(output_stack, local_stack, &a->min, &b->min),
        rational_add(output_stack, local_stack, &a->max, &b->max) };
}

struct RationalInterval rational_interval_negate(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a)
{
    struct RationalInterval out;
    if (rational_get_sign(&a->min) != -rational_get_sign(&a->max))
    {
        out.min = rational_negate(output_stack, &a->max);
        out.max = rational_negate(output_stack, &a->min);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational max_negation = rational_negate(local_stack, &a->max);
        out.min = rational_get_min(output_stack, local_stack, &max_negation, &a->min);
        struct Rational min_negation = rational_negate(local_stack, &a->min);
        out.max = rational_get_max(output_stack, local_stack, &min_negation, &a->max);
        local_stack->cursor = local_stack_savepoint;
    }
    return out;
}

struct RationalInterval rational_interval_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a, struct RationalInterval*b)
{
    struct RationalInterval out;
    if (rational_get_sign(&a->min) >= 0)
    {
        if (rational_get_sign(&b->min) >= 0)
        {
            out.min = rational_multiply(output_stack, local_stack, &a->min, &b->min);
            out.max = rational_multiply(output_stack, local_stack, &a->max, &b->max);
        }
        else
        {
            out.min = rational_multiply(output_stack, local_stack, &a->max, &b->min);
            if (rational_get_sign(&b->max) <= 0)
            {
                out.max = rational_multiply(output_stack, local_stack, &a->min, &b->max);
            }
            else
            {
                out.max = rational_multiply(output_stack, local_stack, &a->max, &b->max);
            }
        }
    }
    else if (rational_get_sign(&a->max) <= 0)
    {
        if (rational_get_sign(&b->min) >= 0)
        {
            out.min = rational_multiply(output_stack, local_stack, &a->min, &b->max);
            out.max = rational_multiply(output_stack, local_stack, &a->max, &b->min);
        }
        else
        {
            if (rational_get_sign(&b->max) <= 0)
            {
                out.min = rational_multiply(output_stack, local_stack, &a->max, &b->max);
            }
            else
            {
                out.min = rational_multiply(output_stack, local_stack, &a->min, &b->max);
            }
            out.max = rational_multiply(output_stack, local_stack, &a->min, &b->min);
        }
    }
    else
    {
        if (rational_get_sign(&b->min) >= 0)
        {
            out.min = rational_multiply(output_stack, local_stack, &a->min, &b->max);
            out.max = rational_multiply(output_stack, local_stack, &a->max, &b->max);
        }
        else if (rational_get_sign(&b->max) <= 0)
        {
            out.min = rational_multiply(output_stack, local_stack, &a->max, &b->min);
            out.max = rational_multiply(output_stack, local_stack, &a->min, &b->min);
        }
        else
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Rational product_a =
                rational_multiply(local_stack, output_stack, &a->min, &b->max);
            struct Rational product_b =
                rational_multiply(local_stack, output_stack, &a->max, &b->min);
            out.min = rational_get_min(output_stack, local_stack, &product_a, &product_b);
            product_a = rational_multiply(local_stack, output_stack, &a->min, &b->min);
            product_b = rational_multiply(local_stack, output_stack, &a->max, &b->max);
            out.max = rational_get_max(output_stack, local_stack, &product_a, &product_b);
            local_stack->cursor = local_stack_savepoint;
        }
    }
    return out;
}

void rational_interval_expand_bounds(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalInterval*a, struct Rational*bound_candidate)
{
    if (rational_compare(local_stack_a, local_stack_b, bound_candidate, &a->min) < 0)
    {
        a->min = *bound_candidate;
    }
    else if (rational_compare(local_stack_a, local_stack_b, bound_candidate, &a->max) > 0)
    {
        a->max = *bound_candidate;
    }
}

struct Rational rational_interval_get_factor_interval_size(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*factor_a,
    struct RationalInterval*factor_b, struct Rational*product_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational a_max_magnitude =
        rational_interval_get_max_magnitude(local_stack, output_stack, factor_a);
    struct Rational b_max_magnitude =
        rational_interval_get_max_magnitude(local_stack, output_stack, factor_b);
    struct Rational max_magnitude_sum =
        rational_add(local_stack, output_stack, &a_max_magnitude, &b_max_magnitude);
    struct Rational divisor =
        rational_integer_add(local_stack, output_stack, &max_magnitude_sum, &one);
    struct Rational out =
        rational_divide(output_stack, local_stack, product_interval_size, &divisor);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Requires that if x straddles an axis, y does not and vice versa.
struct RationalInterval rational_interval_estimate_atan2(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*y, struct RationalInterval*x,
    struct Rational*bound_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out;
    struct RationalInterval bound_estimate;
    if (x->min.numerator->sign >= 0)
    {
        if (y->min.numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->max,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->min,
                *bound_interval_size);
        }
        else if (y->max.numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->min,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->min,
                *bound_interval_size);
        }
        else
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->min,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->max,
                *bound_interval_size);
        }
    }
    else if (x->max.numerator->sign >= 0)
    {
        ASSERT(y->min.numerator->sign == -y->max.numerator->sign,
            "rational_interval_estimate_atan2 x and y arguments both straddled an axis.");
        if (y->min.numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->max,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->min,
                *bound_interval_size);
        }
        else
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->min,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->max,
                *bound_interval_size);
        }
    }
    else
    {
        if (y->min.numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->max,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->min,
                *bound_interval_size);
        }
        else if (y->max.numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->max,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->max,
                *bound_interval_size);
        }
        else
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->max, &x->min,
                *bound_interval_size);
            out.min = rational_copy(output_stack, &bound_estimate.min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, &y->min, &x->max,
                *bound_interval_size);
        }
    }
    out.max = rational_copy(output_stack, &bound_estimate.max);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval rational_interval_to_float_interval(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*a,
    struct Rational*bound_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval out =
    { rational_get_float_estimate(local_stack, output_stack, &a->min, bound_interval_size).min,
    rational_get_float_estimate(local_stack, output_stack, &a->max, bound_interval_size).max };
    out = float_interval_copy(output_stack, &out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

#include "integer/rational/rational_interval/pi/pi.c"