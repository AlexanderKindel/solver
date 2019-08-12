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

struct FloatInterval*float_interval_copy(struct Stack*output_stack, struct FloatInterval*a)
{
    struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
    out->min = float_copy(output_stack, a->min);
    out->max = float_copy(output_stack, a->max);
    return out;
}

bool float_intervals_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval*a, struct FloatInterval*b)
{
    return float_compare(stack_a, stack_b, a->min, b->max) > 0 ||
        float_compare(stack_a, stack_b, b->min, a->max) > 0;
}

struct FloatInterval*float_interval_add(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b)
{
    struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
    out->min = float_add(output_stack, local_stack, a->min, a->min);
    out->max = float_add(output_stack, local_stack, a->max, a->max);
    return out;
}

struct FloatInterval*float_interval_negative(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a)
{
    struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
    if (a->min->significand->sign != -a->max->significand->sign)
    {
        out->min = float_negative(output_stack, a->max);
        out->max = float_negative(output_stack, a->min);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        out->min = float_copy(output_stack, float_min(output_stack, local_stack, a->min,
            float_negative(local_stack, a->max)));
        out->max = float_copy(output_stack, float_max(output_stack, local_stack,
            float_negative(local_stack, a->min), a->max));
        local_stack->cursor = local_stack_savepoint;
    }
    return out;
}

struct FloatInterval*float_interval_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval*out = float_interval_add(output_stack, local_stack, a,
        float_interval_negative(local_stack, output_stack, b));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*float_interval_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a, struct FloatInterval*b)
{
    struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
    if (a->min->significand->sign >= 0)
    {
        if (b->min->significand->sign >= 0)
        {
            out->min = float_multiply(output_stack, local_stack, a->min, b->min);
        }
        else
        {
            out->min = float_multiply(output_stack, local_stack, a->max, b->min);
        }
        out->max = float_multiply(output_stack, local_stack, a->max, b->max);
    }
    else
    {
        if (b->min->significand->sign >= 0)
        {
            out->min = float_multiply(output_stack, local_stack, a->min, b->max);
            out->max = float_multiply(output_stack, local_stack, a->max, b->max);
        }
        else
        {
            void*local_stack_savepoint = local_stack->cursor;
            out->min = float_copy(output_stack, float_min(output_stack, local_stack,
                float_multiply(local_stack, output_stack, a->min, b->max),
                float_multiply(local_stack, output_stack, a->max, b->min)));
            out->max = float_copy(output_stack, float_max(output_stack, local_stack,
                float_multiply(local_stack, output_stack, a->min, b->min),
                float_multiply(local_stack, output_stack, a->max, b->max)));
            local_stack->cursor = local_stack_savepoint;
        }
    }
    return out;
}

void float_interval_to_rational_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct FloatInterval*a)
{
    out->min = float_to_rational(output_stack, local_stack, a->min);
    out->max = float_to_rational(output_stack, local_stack, a->max);
}

bool rectangular_estimates_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct RectangularEstimate*a, struct RectangularEstimate*b)
{
    return float_intervals_are_disjoint(stack_a, stack_b, a->real_part_estimate,
        b->real_part_estimate) ||
        float_intervals_are_disjoint(stack_a, stack_b, a->imaginary_part_estimate,
            b->imaginary_part_estimate);
}