#include "declarations.h"

struct Rational*rational_interval_max_magnitude(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_max(output_stack, local_stack,
        rational_magnitude(local_stack, a->min), rational_magnitude(local_stack, a->max));
    out = rational_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void rational_interval_expand_bounds(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalInterval*a, struct Rational*bound_candidate)
{
    if (rational_compare(stack_a, stack_b, bound_candidate, a->min) < 0)
    {
        a->min = bound_candidate;
    }
    else if (rational_compare(stack_a, stack_b, bound_candidate, a->max) > 0)
    {
        a->max = bound_candidate;
    }
}

struct Rational*rational_interval_factor_interval_size(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*factor_a, struct RationalInterval*factor_b,
    struct Rational*product_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_divide(output_stack, local_stack, product_interval_size,
        rational_add(local_stack, output_stack,
            rational_interval_max_magnitude(local_stack, output_stack, factor_a),
            rational_integer_add(local_stack, output_stack,
                rational_interval_max_magnitude(local_stack, output_stack, factor_b), &one)));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Requires that if x straddles an axis, y does not and vice versa.
void rational_interval_estimate_atan2(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct RationalInterval*y, struct RationalInterval*x,
    struct Rational*bound_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval bound_estimate;
    if (x->min->numerator->sign >= 0)
    {
        if (y->min->numerator->sign >= 0)
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->max,
                bound_interval_size);
        }
        else if (y->max->numerator->sign >= 0)
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->max,
                bound_interval_size);
        }
        else
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->max,
                bound_interval_size);
        }
    }
    else if (x->max->numerator->sign >= 0)
    {
        ASSERT(y->min->numerator->sign == -y->max->numerator->sign,
            "rational_interval_estimate_atan2 x and y arguments both straddled an axis.");
        if (y->min->numerator->sign >= 0)
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->min,
                bound_interval_size);
        }
        else
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->max,
                bound_interval_size);
        }
    }
    else
    {
        if (y->min->numerator->sign >= 0)
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->min,
                bound_interval_size);
        }
        else if (y->max->numerator->sign >= 0)
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->min,
                bound_interval_size);
        }
        else
        {
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->min, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate.min);
            rational_estimate_atan2(local_stack, output_stack, &bound_estimate, x->max, y->min,
                bound_interval_size);
        }
    }
    out->max = rational_copy(output_stack, bound_estimate.max);
    local_stack->cursor = local_stack_savepoint;
}

void rational_interval_to_float_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct RationalInterval*a, struct Rational*bound_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval min_estimate;
    rational_float_estimate(local_stack, output_stack, &min_estimate, a->min, bound_interval_size);
    rational_float_estimate(local_stack, output_stack, out, a->max, bound_interval_size);
    out->min = float_copy(output_stack, min_estimate.min);
    out->max = float_copy(output_stack, out->max);
    local_stack->cursor = local_stack_savepoint;
}

bool float_intervals_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval*a, struct FloatInterval*b)
{
    return float_compare(stack_a, stack_b, a->min, b->max) > 0 ||
        float_compare(stack_a, stack_b, b->min, a->max) > 0;
}

void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct FloatInterval*a)
{
    out->min = float_to_rational(output_stack, local_stack, a->min);
    out->max = float_to_rational(output_stack, local_stack, a->max);
}