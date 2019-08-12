#include "declarations.h"

struct Rational*rational_copy(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    out->numerator = integer_copy(output_stack, a->numerator);
    out->denominator = integer_copy(output_stack, a->denominator);
    return out;
}

struct Rational*rational_reduced(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*numerator, struct Integer*denominator)
{
    if (denominator->sign == 0)
    {
        puts("Tried to divide by 0.");
        return 0;
    }
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    if (numerator->sign == 0)
    {
        out->numerator = &zero;
        out->denominator = &one;
        return out;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo gcd_info;
    integer_extended_gcd(local_stack, output_stack, &gcd_info, numerator, denominator);
    out->numerator = gcd_info.a_over_gcd;
    out->denominator = gcd_info.b_over_gcd;
    out->numerator->sign *= out->denominator->sign;
    out->denominator->sign *= out->denominator->sign;
    out->numerator = integer_copy(output_stack, out->numerator);
    out->denominator = integer_copy(output_stack, out->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

bool rational_equals(struct Rational*a, struct Rational*b)
{
    return integer_equals(a->numerator, b->numerator) &&
        integer_equals(a->denominator, b->denominator);
}

struct Rational*rational_magnitude(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    out->numerator = integer_magnitude(output_stack, a->numerator);
    out->denominator = integer_copy(output_stack, a->denominator);
    return out;
}

struct Rational*rational_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack, integer_add(local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b->denominator),
        integer_multiply(local_stack, output_stack, b->numerator, a->denominator)),
        integer_multiply(local_stack, output_stack, a->denominator, b->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_integer_add(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    out->numerator = integer_add(output_stack, a->numerator,
        integer_multiply(local_stack, output_stack, b, a->denominator));
    out->denominator = integer_copy(output_stack, a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_negative(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    out->numerator = integer_negative(output_stack, a->numerator);
    out->denominator = integer_copy(output_stack, a->denominator);
    return out;
}

struct Rational*rational_subtract(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*minuend, struct Rational*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_subtract(local_stack, output_stack,
            integer_multiply(local_stack, output_stack, minuend->numerator,
                subtrahend->denominator),
            integer_multiply(local_stack, output_stack, subtrahend->numerator,
                minuend->denominator)),
        integer_multiply(local_stack, output_stack, minuend->denominator,
            subtrahend->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b->numerator),
        integer_multiply(local_stack, output_stack, a->denominator, b->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_unreduced_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Rational*b, void*unused)
{
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    out->numerator = integer_multiply(output_stack, local_stack, a->numerator, b->numerator);
    out->denominator = integer_multiply(output_stack, local_stack, a->denominator, b->denominator);
    return out;
}

struct Rational*rational_integer_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b), a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_reciprocal(struct Stack*output_stack, struct Rational*a)
{
    struct Rational*out = ALLOCATE(output_stack, struct Rational);
    out->numerator = integer_copy(output_stack, a->denominator);
    out->denominator = integer_copy(output_stack, a->numerator);
    out->numerator->sign *= out->denominator->sign;
    out->denominator->sign = 1;
    return out;
}

struct Rational*rational_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Rational*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, dividend->numerator, divisor->denominator),
        integer_multiply(local_stack, output_stack, dividend->denominator, divisor->numerator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_integer_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*dividend, struct Integer*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack, dividend->numerator,
        integer_multiply(local_stack, output_stack, dividend->denominator, divisor));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_doubled(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack,
        integer_doubled(local_stack, a->numerator), a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_half(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_reduced(output_stack, local_stack, a->numerator,
        integer_half(local_stack, a->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

int8_t rational_compare(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t out = integer_compare(stack_a, stack_b,
        integer_multiply(stack_a, stack_b, a->numerator, b->denominator),
        integer_multiply(stack_a, stack_b, a->denominator, b->numerator));
    stack_a->cursor = stack_a_savepoint;
    return out;
}

struct Rational*rational_min(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    if (rational_compare(stack_a, stack_b, a, b) < 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

struct Rational*rational_max(struct Stack*stack_a, struct Stack*stack_b, struct Rational*a,
    struct Rational*b)
{
    if (rational_compare(stack_a, stack_b, a, b) > 0)
    {
        return a;
    }
    else
    {
        return b;
    }
}

struct Rational*rational_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*base, struct Integer*exponent)
{
    return generic_exponentiate(&(struct RingOperations) { rational_copy, rational_equals,
        &rational_zero, &rational_one, rational_generic_add, rational_generic_negative,
        rational_unreduced_multiply }, output_stack, local_stack, base, exponent, 0);
}

struct RationalPolynomial*rational_minimal_polynomial(struct Stack*output_stack, struct Rational*a)
{
    struct RationalPolynomial*out = polynomial_allocate(output_stack, 2);
    out->coefficients[0] = ALLOCATE(output_stack, struct Rational);
    out->coefficients[0]->numerator = integer_negative(output_stack, a->numerator);
    out->coefficients[0]->denominator = integer_copy(output_stack, a->denominator);
    out->coefficients[1] = &rational_one;
    return out;
}

void rational_estimate_size_or_cosine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a_squared, struct Integer*factorial_component,
    struct Rational*delta, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (rational_compare(output_stack, local_stack, rational_magnitude(local_stack, delta),
        interval_size) > 0)
    {
        out->min = rational_add(local_stack, output_stack, out->min, delta);
        factorial_component = integer_add(local_stack, factorial_component, &one);
        delta = rational_multiply(local_stack, output_stack, delta,
            rational_integer_divide(local_stack, output_stack, a_squared, factorial_component));
        factorial_component = integer_add(local_stack, factorial_component, &one);
        struct Integer*negative_factorial_component =
            integer_negative(local_stack, factorial_component);
        delta =
            rational_integer_divide(local_stack, output_stack, delta, negative_factorial_component);
    }
    out->min = rational_copy(output_stack, out->min);
    if (delta->numerator->value_count > 0)
    {
        out->max = rational_add(output_stack, local_stack, out->min, delta);
    }
    else
    {
        out->max = out->min;
        out->min = rational_add(output_stack, local_stack, out->max, delta);
    }
    local_stack->cursor = local_stack_savepoint;
}

void rational_estimate_cosine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->min->numerator = &one;
    out->min->denominator = &one;
    struct Rational*a_squared = rational_multiply(local_stack, output_stack, a, a);
    struct Integer*factorial_component = &INT(2, +);
    struct Rational*delta =
        rational_integer_divide(local_stack, output_stack, a_squared, &INT(2, -));
    rational_estimate_size_or_cosine(output_stack, local_stack, out, a_squared, factorial_component,
        delta, interval_size);
    local_stack->cursor = local_stack_savepoint;
}

void rational_estimate_sine(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->min = a;
    struct Rational*a_squared = rational_multiply(local_stack, output_stack, a, a);
    struct Integer*factorial_component = &INT(3, +);
    struct Rational*delta = rational_integer_divide(local_stack, output_stack,
        rational_multiply(local_stack, output_stack, a_squared, a), &INT(6, -));
    rational_estimate_size_or_cosine(output_stack, local_stack, out, a_squared, factorial_component,
        delta, interval_size);
    local_stack->cursor = local_stack_savepoint;
}

void rational_estimate_arctangent(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    ASSERT(rational_compare(output_stack, local_stack, rational_magnitude(local_stack, a),
        &rational_one) <= 0, "rational_estimate_arctangent was called with an a argument whose "
        "magnitude was greater than one. Use an identity to evaluate an equivalent expression "
        "involving its reciprocal instead.");
    out->min = a;
    struct Rational*negative_a_squared = rational_multiply(local_stack, output_stack, a, a);
    negative_a_squared->numerator->sign *= -1;
    struct Rational*delta_numerator =
        rational_multiply(local_stack, output_stack, negative_a_squared, a);
    struct Integer*delta_denominator = &INT(3, +);
    while (true)
    {
        struct Rational*delta =
            rational_integer_divide(local_stack, output_stack, delta_numerator, delta_denominator);
        if (rational_compare(output_stack, local_stack, rational_magnitude(local_stack, delta),
            interval_size) <= 0)
        {
            out->min = rational_copy(output_stack, out->min);
            if (delta->numerator->value_count > 0)
            {
                out->max = rational_add(output_stack, local_stack, out->min, delta);
            }
            else
            {
                out->max = out->min;
                out->min = rational_add(output_stack, local_stack, out->max, delta);
            }
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        out->min = rational_add(local_stack, output_stack, out->min, delta);
        delta_numerator =
            rational_multiply(local_stack, output_stack, delta_numerator, negative_a_squared);
        delta_denominator = integer_add(local_stack, delta_denominator, &INT(2, +));
    }
}

void rational_estimate_atan2(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Rational*y, struct Rational*x,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    if (x->numerator->sign == 0)
    {
        struct Integer pi_multiple_numerator;
        if (y->numerator->sign > 0)
        {
            pi_multiple_numerator = one;
        }
        else
        {
            pi_multiple_numerator = INT(3, +);
        }
        struct Rational pi_multiple = { &pi_multiple_numerator, &INT(2, +) };
        interval_size = rational_divide(local_stack, output_stack, interval_size, &pi_multiple);
        pi_estimate(interval_size);
        out->min = rational_multiply(output_stack, local_stack, &pi_multiple, pi.min);
        out->max = rational_multiply(output_stack, local_stack, &pi_multiple, pi.max);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    struct Rational*ratio = rational_divide(local_stack, output_stack, y, x);
    int8_t ratio_magnitude_comparison = rational_compare(output_stack, local_stack,
        rational_magnitude(local_stack, ratio), &rational_one);
    if (ratio_magnitude_comparison == 0)
    {
        struct Integer pi_multiple_numerator;
        if (x->numerator->sign > 0)
        {
            if (y->numerator->sign > 0)
            {
                pi_multiple_numerator = one;
            }
            else
            {
                pi_multiple_numerator = INT(7, +);
            }
        }
        else
        {
            if (y->numerator->sign > 0)
            {
                pi_multiple_numerator = INT(3, +);
            }
            else
            {
                pi_multiple_numerator = INT(5, +);
            }
        }
        struct Rational pi_multiple = { &pi_multiple_numerator, &INT(4, +) };
        interval_size = rational_divide(local_stack, output_stack, interval_size, &pi_multiple);
        pi_estimate(interval_size);
        out->min = rational_multiply(output_stack, local_stack, &pi_multiple, pi.min);
        out->max = rational_multiply(output_stack, local_stack, &pi_multiple, pi.max);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    struct Rational*multiple_of_pi_to_add = &(struct Rational) { &INT(2, +), &one };
    if (x->numerator->sign < 0)
    {
        multiple_of_pi_to_add = &rational_one;
    }
    else if (y->numerator->sign >= 0)
    {
        multiple_of_pi_to_add = &rational_zero;
    }
    interval_size = rational_half(local_stack, output_stack, interval_size);
    if (ratio_magnitude_comparison < 0)
    {
        rational_estimate_arctangent(local_stack, output_stack, out, ratio, interval_size);
    }
    else
    {
        if (ratio->numerator->sign > 0)
        {
            multiple_of_pi_to_add = rational_add(local_stack, output_stack, multiple_of_pi_to_add,
                &(struct Rational){&one, &INT(2, +)});
        }
        else
        {
            multiple_of_pi_to_add = rational_add(local_stack, output_stack, multiple_of_pi_to_add,
                &(struct Rational){&INT(1, -), &INT(2, +)});
        }
        struct RationalInterval arctan_of_reciprocal;
        rational_estimate_arctangent(local_stack, output_stack, &arctan_of_reciprocal,
            rational_reciprocal(local_stack, ratio), interval_size);
        out->min = arctan_of_reciprocal.max;
        out->min->numerator->sign *= -1;
        out->max = arctan_of_reciprocal.min;
        out->max->numerator->sign *= -1;
    }
    if (multiple_of_pi_to_add->numerator->sign == 0)
    {
        out->min = rational_copy(output_stack, out->min);
        out->max = rational_copy(output_stack, out->max);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    pi_estimate(rational_divide(local_stack, output_stack, interval_size,
        rational_magnitude(local_stack, multiple_of_pi_to_add)));
    struct Rational*pi_estimate_min_multiple =
        rational_multiply(local_stack, output_stack, multiple_of_pi_to_add, pi.min);
    struct Rational*pi_estimate_max_multiple =
        rational_multiply(local_stack, output_stack, multiple_of_pi_to_add, pi.max);
    if (multiple_of_pi_to_add->numerator->sign > 0)
    {
        out->min = rational_add(output_stack, local_stack, out->min, pi_estimate_min_multiple);
        out->max = rational_add(output_stack, local_stack, out->max, pi_estimate_max_multiple);
    }
    else
    {
        out->min = rational_add(output_stack, local_stack, out->min, pi_estimate_max_multiple);
        out->max = rational_add(output_stack, local_stack, out->max, pi_estimate_min_multiple);
    }
    local_stack->cursor = local_stack_savepoint;
}

void positive_rational_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerDivision division;
    integer_euclidean_divide(local_stack, output_stack, &division, a->numerator, a->denominator);
    (*out_min)->significand = division.quotient;
    (*out_min)->exponent = &zero;
    struct Integer*estimate_denominator = &one;
    while (true)
    {
        if (division.remainder->value_count == 0)
        {
            *out_min = float_copy(output_stack, *out_min);
            *out_max = *out_min;
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        if (integer_compare(output_stack, local_stack, interval_size->denominator,
            integer_multiply(local_stack, output_stack, estimate_denominator,
                interval_size->numerator)) <= 0)
        {
            break;
        }
        integer_euclidean_divide(local_stack, output_stack, &division,
            integer_doubled(local_stack, division.remainder), a->denominator);
        estimate_denominator = integer_doubled(local_stack, estimate_denominator);
        (*out_min)->significand = integer_add(local_stack, division.quotient,
            integer_doubled(local_stack, (*out_min)->significand));
        (*out_min)->exponent = integer_add(local_stack, (*out_min)->exponent, &INT(1, -));
    }
    *out_min =
        float_reduced(output_stack, local_stack, (*out_min)->significand, (*out_min)->exponent);
    *out_max = float_reduced(output_stack, local_stack,
        integer_add(local_stack, (*out_min)->significand, &one), (*out_min)->exponent);
    local_stack->cursor = local_stack_savepoint;
}

void rational_float_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->min = ALLOCATE(output_stack, struct Float);
    out->max = ALLOCATE(output_stack, struct Float);
    if (a->numerator->sign < 0)
    {
        positive_rational_float_estimate(output_stack, local_stack, &out->max, &out->min,
            rational_negative(local_stack, a), interval_size);
        out->min->significand->sign = -1;
        out->max->significand->sign = -1;
    }
    else
    {
        positive_rational_float_estimate(output_stack, local_stack, &out->min, &out->max, a,
            interval_size);
    }
    local_stack->cursor = local_stack_savepoint;
}