#include "declarations.h"

void*interval_max_magnitude(struct IntervalBoundOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Interval*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    void*out = operations->bmax(output_stack, local_stack,
        operations->magnitude(local_stack, a->min), operations->magnitude(local_stack, a->max));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*interval_add(struct IntervalBoundOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Interval*a, struct Interval*b)
{
    struct Interval*out = ALLOCATE(output_stack, struct Interval);
    out->min = operations->add(output_stack, local_stack, a->min, b->min);
    out->max = operations->add(output_stack, local_stack, a->max, b->max);
    return out;
}

void*interval_negative(struct IntervalBoundOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Interval*a)
{
    struct Interval*out = ALLOCATE(output_stack, struct Interval);
    if (operations->sign(a->min) != -operations->sign(a->max))
    {
        out->min = operations->negative(output_stack, a->max);
        out->max = operations->negative(output_stack, a->min);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        out->min = operations->copy(output_stack, operations->bmin(output_stack, local_stack,
            operations->negative(local_stack, a->max), a->min));
        out->max = operations->copy(output_stack, operations->bmax(output_stack, local_stack,
            operations->negative(local_stack, a->min), a->max));
        local_stack->cursor = local_stack_savepoint;
    }
    return out;
}

void*interval_multiply(struct IntervalBoundOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Interval*a, struct Interval*b)
{
    struct Interval*out = ALLOCATE(output_stack, struct Interval);
    if (operations->sign(a->min) >= 0)
    {
        if (operations->sign(b->min) >= 0)
        {
            out->min = operations->multiply(output_stack, local_stack, a->min, b->min);
            out->max = operations->multiply(output_stack, local_stack, a->max, b->max);
        }
        else
        {
            out->min = operations->multiply(output_stack, local_stack, a->max, b->min);
            if (operations->sign(b->max) <= 0)
            {
                out->max = operations->multiply(output_stack, local_stack, a->min, b->max);
            }
            else
            {
                out->max = operations->multiply(output_stack, local_stack, a->max, b->max);
            }
        }
    }
    else if (operations->sign(a->max) <= 0)
    {
        if (operations->sign(b->min) >= 0)
        {
            out->min = operations->multiply(output_stack, local_stack, a->min, b->max);
            out->max = operations->multiply(output_stack, local_stack, a->max, b->min);
        }
        else
        {
            if (operations->sign(b->max) <= 0)
            {
                out->min = operations->multiply(output_stack, local_stack, a->max, b->max);
            }
            else
            {
                out->min = operations->multiply(output_stack, local_stack, a->min, b->max);
            }
            out->max = operations->multiply(output_stack, local_stack, a->min, b->min);
        }
    }
    else
    {
        if (operations->sign(b->min) >= 0)
        {
            out->min = operations->multiply(output_stack, local_stack, a->min, b->max);
            out->max = operations->multiply(output_stack, local_stack, a->max, b->max);
        }
        else if (operations->sign(b->max) <= 0)
        {
            out->min = operations->multiply(output_stack, local_stack, a->max, b->min);
            out->max = operations->multiply(output_stack, local_stack, a->min, b->min);
        }
        else
        {
            void*local_stack_savepoint = local_stack->cursor;
            out->min = operations->copy(output_stack, operations->bmin(output_stack, local_stack,
                operations->multiply(local_stack, output_stack, a->min, b->max),
                operations->multiply(local_stack, output_stack, a->max, b->min)));
            out->max = operations->copy(output_stack, operations->bmax(output_stack, local_stack,
                operations->multiply(local_stack, output_stack, a->min, b->min),
                operations->multiply(local_stack, output_stack, a->max, b->max)));
            local_stack->cursor = local_stack_savepoint;
        }
    }
    return out;
}

struct RationalInterval*rational_interval_copy(struct Stack*output_stack, struct RationalInterval*a)
{
    struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
    out->min = rational_copy(output_stack, a->min);
    out->max = rational_copy(output_stack, a->max);
    return out;
}

struct Rational*rational_interval_max_magnitude(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*a)
{
    return interval_max_magnitude(&rational_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a);
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

struct RationalInterval*rational_interval_add(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a, struct RationalInterval*b)
{
    return interval_add(&rational_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a, (struct Interval*)b);
}

struct RationalInterval*rational_interval_negative(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a)
{
    return interval_negative(&rational_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a);
}

struct RationalInterval*rational_interval_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a, struct RationalInterval*b)
{
    return interval_multiply(&rational_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a, (struct Interval*)b);
}

//Requires that if x straddles an axis, y does not and vice versa.
struct RationalInterval*rational_interval_estimate_atan2(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*y, struct RationalInterval*x,
    struct Rational*bound_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
    struct RationalInterval*bound_estimate;
    if (x->min->numerator->sign >= 0)
    {
        if (y->min->numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->max,
                bound_interval_size);
        }
        else if (y->max->numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->max,
                bound_interval_size);
        }
        else
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->max,
                bound_interval_size);
        }
    }
    else if (x->max->numerator->sign >= 0)
    {
        ASSERT(y->min->numerator->sign == -y->max->numerator->sign,
            "rational_interval_estimate_atan2 x and y arguments both straddled an axis.");
        if (y->min->numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->min,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->min,
                bound_interval_size);
        }
        else
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->max,
                bound_interval_size);
        }
    }
    else
    {
        if (y->min->numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->min,
                bound_interval_size);
        }
        else if (y->max->numerator->sign >= 0)
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->min,
                bound_interval_size);
        }
        else
        {
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->min, y->max,
                bound_interval_size);
            out->min = rational_copy(output_stack, bound_estimate->min);
            bound_estimate = rational_estimate_atan2(local_stack, output_stack, x->max, y->min,
                bound_interval_size);
        }
    }
    out->max = rational_copy(output_stack, bound_estimate->max);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*rational_interval_to_float_interval(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*a, struct Rational*bound_interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
    out->min = float_copy(output_stack,
        rational_float_estimate(local_stack, output_stack, a->min, bound_interval_size)->min);
    out->max = float_copy(output_stack,
        rational_float_estimate(local_stack, output_stack, a->max, bound_interval_size)->max);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval*float_interval_copy(struct Stack*output_stack, struct FloatInterval*a)
{
    struct FloatInterval*out = ALLOCATE(output_stack, struct FloatInterval);
    out->min = float_copy(output_stack, a->min);
    out->max = float_copy(output_stack, a->max);
    return out;
}

struct Float*float_interval_max_magnitude(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a)
{
    return interval_max_magnitude(&float_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a);
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
    return interval_add(&float_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a, (struct Interval*)b);
}

struct FloatInterval*float_interval_negative(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*a)
{
    return interval_negative(&float_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a);
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
    return interval_multiply(&float_interval_bound_operations, output_stack, local_stack,
        (struct Interval*)a, (struct Interval*)b);
}

struct RationalInterval*float_interval_to_rational_interval(struct Stack*output_stack,
    struct Stack*local_stack, struct FloatInterval*a)
{
    struct RationalInterval*out = ALLOCATE(output_stack, struct RationalInterval);
    out->min = float_to_rational(output_stack, local_stack, a->min);
    out->max = float_to_rational(output_stack, local_stack, a->max);
    return out;
}

bool rectangular_estimates_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct RectangularEstimate*a, struct RectangularEstimate*b)
{
    return float_intervals_are_disjoint(stack_a, stack_b, a->real_part_estimate,
        b->real_part_estimate) ||
        float_intervals_are_disjoint(stack_a, stack_b, a->imaginary_part_estimate,
            b->imaginary_part_estimate);
}