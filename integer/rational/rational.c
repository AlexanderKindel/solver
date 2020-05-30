#include "declarations.h"

struct Rational rational_copy(struct Stack*output_stack, struct Rational*a)
{
    return (struct Rational) { integer_copy(output_stack, a->numerator),
        integer_copy(output_stack, a->denominator) };
}

struct Rational rational_reduce(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Integer*numerator, struct Integer*denominator)
{
    ASSERT(denominator->value_count, "rational_reduced denominator was 0.");
    if (numerator->sign == 0)
    {
        return (struct Rational) { &zero, &one };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo gcd_info =
        integer_get_extended_gcd(local_stack, output_stack, numerator, denominator);
    struct Rational out = { integer_copy(output_stack, gcd_info.a_over_gcd),
        integer_copy(output_stack, gcd_info.b_over_gcd) };
    out.numerator->sign *= out.denominator->sign;
    out.denominator->sign *= out.denominator->sign;
    local_stack->cursor = local_stack_savepoint;
    return out;
}

bool rational_equals(struct Rational*a, struct Rational*b)
{
    return integer_equals(a->numerator, b->numerator) &&
        integer_equals(a->denominator, b->denominator);
}

struct Rational rational_get_magnitude(struct Stack*output_stack, struct Rational*a)
{
    return (struct Rational) { integer_get_magnitude(output_stack, a->numerator),
        integer_copy(output_stack, a->denominator) };
}

struct Rational rational_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Rational*a, struct Rational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack, integer_add(local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b->denominator),
        integer_multiply(local_stack, output_stack, b->numerator, a->denominator)),
        integer_multiply(local_stack, output_stack, a->denominator, b->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_integer_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = { integer_add(output_stack, a->numerator,
        integer_multiply(local_stack, output_stack, b, a->denominator)),
        integer_copy(output_stack, a->denominator) };
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_negate(struct Stack*output_stack, struct Rational*a)
{
    return (struct Rational) { integer_negate(output_stack, a->numerator),
        integer_copy(output_stack, a->denominator) };
}

struct Rational rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*minuend, struct Rational*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack,
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

struct Rational rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b->numerator),
        integer_multiply(local_stack, output_stack, a->denominator, b->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_unreduced_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b)
{
    return (struct Rational) {
        integer_multiply(output_stack, local_stack, a->numerator, b->numerator),
        integer_multiply(output_stack, local_stack, a->denominator, b->denominator) };
}

struct Rational rational_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Integer*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, a->numerator, b), a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_get_reciprocal(struct Stack*output_stack, struct Rational*a)
{
    struct Rational out =
    { integer_copy(output_stack, a->denominator), integer_copy(output_stack, a->numerator) };
    out.numerator->sign *= out.denominator->sign;
    out.denominator->sign = 1;
    return out;
}

struct Rational rational_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*dividend, struct Rational*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack,
        integer_multiply(local_stack, output_stack, dividend->numerator, divisor->denominator),
        integer_multiply(local_stack, output_stack, dividend->denominator, divisor->numerator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_integer_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*dividend, struct Integer*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack, dividend->numerator,
        integer_multiply(local_stack, output_stack, dividend->denominator, divisor));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_double(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack,
        integer_double(local_stack, a->numerator), a->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_halve(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_reduce(output_stack, local_stack, a->numerator,
        integer_double(local_stack, a->denominator));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational rational_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational base, struct Integer*exponent)
{
    if (!exponent->value_count)
    {
        return rational_one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_one;
    while (true)
    {
        if (exponent->value[0] & 1)
        {
            out = rational_unreduced_multiply(local_stack, output_stack, &out, &base);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = rational_copy(output_stack, &out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = rational_unreduced_multiply(local_stack, output_stack, &base, &base);
    }
}

int8_t rational_get_sign(struct Rational*a)
{
    return a->numerator->sign;
}

int8_t rational_compare(struct Stack*restrict local_stack_a, struct Stack*restrict local_stack_b,
    struct Rational*a, struct Rational*b)
{
    void*local_stack_a_savepoint = local_stack_a->cursor;
    int8_t out = integer_compare(local_stack_a, local_stack_b,
        integer_multiply(local_stack_a, local_stack_b, a->numerator, b->denominator),
        integer_multiply(local_stack_a, local_stack_b, a->denominator, b->numerator));
    local_stack_a->cursor = local_stack_a_savepoint;
    return out;
}

struct Rational rational_get_min(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b)
{
    if (rational_compare(output_stack, local_stack, a, b) < 0)
    {
        return rational_copy(output_stack, a);
    }
    else
    {
        return rational_copy(output_stack, b);
    }
}

struct Rational rational_get_max(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*b)
{
    if (rational_compare(output_stack, local_stack, a, b) > 0)
    {
        return rational_copy(output_stack, a);
    }
    else
    {
        return rational_copy(output_stack, b);
    }
}

//In revolutions.
struct Rational rational_get_argument(struct Stack*output_stack, struct Rational*a)
{
    if (a->numerator->sign < 0)
    {
        return (struct Rational) { &one, integer_initialize(output_stack, 2, 1) };
    }
    return rational_zero;
}

void rational_estimate_sine_or_cosine(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalInterval*out, struct Rational*a_squared,
    struct Integer*factorial_component, struct Rational*delta, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (true)
    {
        struct Rational delta_magnitude = rational_get_magnitude(local_stack, delta);
        if (rational_compare(output_stack, local_stack, &delta_magnitude, interval_size) <= 0)
        {
            break;
        }
        out->min = rational_add(local_stack, output_stack, &out->min, delta);
        factorial_component = integer_add(local_stack, factorial_component, &one);
        struct Rational quotient =
            rational_integer_divide(local_stack, output_stack, a_squared, factorial_component);
        *delta = rational_multiply(local_stack, output_stack, delta, &quotient);
        factorial_component = integer_add(local_stack, factorial_component, &one);
        struct Integer*negative_factorial_component =
            integer_negate(local_stack, factorial_component);
        *delta =
            rational_integer_divide(local_stack, output_stack, delta, negative_factorial_component);
    }
    out->min = rational_copy(output_stack, &out->min);
    if (delta->numerator->value_count > 0)
    {
        out->max = rational_add(output_stack, local_stack, &out->min, delta);
    }
    else
    {
        out->max = out->min;
        out->min = rational_add(output_stack, local_stack, &out->max, delta);
    }
    local_stack->cursor = local_stack_savepoint;
}

struct RationalInterval rational_estimate_cosine(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out;
    out.min = rational_one;
    struct Rational a_squared = rational_multiply(local_stack, output_stack, a, a);
    struct Integer*factorial_component = INT(2, 1);
    struct Rational delta =
        rational_integer_divide(local_stack, output_stack, &a_squared, INT(2, -1));
    rational_estimate_sine_or_cosine(output_stack, local_stack, &out, &a_squared,
        factorial_component, &delta, interval_size);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalInterval rational_estimate_sine(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out;
    out.min = *a;
    struct Rational a_squared = rational_multiply(local_stack, output_stack, a, a);
    struct Integer*factorial_component = INT(3, 1);
    struct Rational product = rational_multiply(local_stack, output_stack, &a_squared, a);
    struct Rational delta =
        rational_integer_divide(local_stack, output_stack, &product, INT(6, -1));
    rational_estimate_sine_or_cosine(output_stack, local_stack, &out, &a_squared,
        factorial_component, &delta, interval_size);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Invalid if the magnitude of the argument is greater than one. In those cases, use an identity to
//evaluate an equivalent expression involving the desired argument's reciprocal instead.
struct RationalInterval rational_estimate_arctangent(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval out;
    out.min = *a;
    struct Rational negative_a_squared = rational_multiply(local_stack, output_stack, a, a);
    negative_a_squared.numerator->sign *= -1;
    struct Rational delta_numerator =
        rational_multiply(local_stack, output_stack, &negative_a_squared, a);
    struct Integer*delta_denominator = INT(3, 1);
    while (true)
    {
        struct Rational delta =
            rational_integer_divide(local_stack, output_stack, &delta_numerator, delta_denominator);
        struct Rational delta_magnitude = rational_get_magnitude(local_stack, &delta);
        if (rational_compare(output_stack, local_stack, &delta_magnitude, interval_size) <= 0)
        {
            out.min = rational_copy(output_stack, &out.min);
            if (delta.numerator->value_count > 0)
            {
                out.max = rational_add(output_stack, local_stack, &out.min, &delta);
            }
            else
            {
                out.max = out.min;
                out.min = rational_add(output_stack, local_stack, &out.max, &delta);
            }
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        out.min = rational_add(local_stack, output_stack, &out.min, &delta);
        delta_numerator =
            rational_multiply(local_stack, output_stack, &delta_numerator, &negative_a_squared);
        delta_denominator = integer_add(local_stack, delta_denominator, INT(2, 1));
    }
}

struct RationalInterval rational_estimate_atan2(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*y, struct Rational*x,
    struct Rational interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    if (x->numerator->sign == 0)
    {
        struct Rational pi_multiple;
        if (y->numerator->sign > 0)
        {
            pi_multiple.numerator = &one;
        }
        else
        {
            pi_multiple.numerator = INT(3, 1);
        }
        pi_multiple.denominator = INT(2, 1);
        interval_size = rational_divide(local_stack, output_stack, &interval_size, &pi_multiple);
        pi_estimate(&interval_size);
        struct RationalInterval out =
        { rational_multiply(output_stack, local_stack, &pi_multiple, &pi.min),
            rational_multiply(output_stack, local_stack, &pi_multiple, &pi.max) };
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Rational ratio = rational_divide(local_stack, output_stack, y, x);
    struct Rational ratio_magnitude = rational_get_magnitude(local_stack, &ratio);
    int8_t ratio_magnitude_comparison =
        rational_compare(output_stack, local_stack, &ratio_magnitude, &rational_one);
    if (ratio_magnitude_comparison == 0)
    {
        struct Rational pi_multiple;
        if (x->numerator->sign > 0)
        {
            if (y->numerator->sign > 0)
            {
                pi_multiple.numerator = &one;
            }
            else
            {
                pi_multiple.numerator = INT(7, 1);
            }
        }
        else
        {
            if (y->numerator->sign > 0)
            {
                pi_multiple.numerator = INT(3, 1);
            }
            else
            {
                pi_multiple.numerator = INT(5, 1);
            }
        }
        pi_multiple.denominator = INT(4, 1);
        interval_size = rational_divide(local_stack, output_stack, &interval_size, &pi_multiple);
        pi_estimate(&interval_size);
        struct RationalInterval out =
        { rational_multiply(output_stack, local_stack, &pi_multiple, &pi.min),
            rational_multiply(output_stack, local_stack, &pi_multiple, &pi.max) };
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Rational multiple_of_pi_to_add = { INT(2, 1), &one };
    if (x->numerator->sign < 0)
    {
        multiple_of_pi_to_add = rational_one;
    }
    else if (y->numerator->sign >= 0)
    {
        multiple_of_pi_to_add = rational_zero;
    }
    interval_size = rational_halve(local_stack, output_stack, &interval_size);
    struct RationalInterval out;
    if (ratio_magnitude_comparison < 0)
    {
        out = rational_estimate_arctangent(local_stack, output_stack, &ratio, &interval_size);
    }
    else
    {
        if (ratio.numerator->sign > 0)
        {
            multiple_of_pi_to_add = rational_add(local_stack, output_stack, &multiple_of_pi_to_add,
                &(struct Rational) { &one, INT(2, 1) });
        }
        else
        {
            multiple_of_pi_to_add = rational_add(local_stack, output_stack, &multiple_of_pi_to_add,
                &(struct Rational) { INT(1, -1), INT(2, 1) });
        }
        struct Rational ratio_reciprocal = rational_get_reciprocal(local_stack, &ratio);
        struct RationalInterval reciprocal_arctangent = rational_estimate_arctangent(local_stack,
            output_stack, &ratio_reciprocal, &interval_size);
        out = rational_interval_negate(local_stack, output_stack, &reciprocal_arctangent);
    }
    if (multiple_of_pi_to_add.numerator->sign == 0)
    {
        out = rational_interval_copy(output_stack, &out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Rational multiple_magnitude =
        rational_get_magnitude(local_stack, &multiple_of_pi_to_add);
    struct Rational magnitude_quotient =
        rational_divide(local_stack, output_stack, &interval_size, &multiple_magnitude);
    pi_estimate(&magnitude_quotient);
    struct RationalInterval pi_estimate_multiple =
    { rational_multiply(local_stack, output_stack, &multiple_of_pi_to_add, &pi.min),
    rational_multiply(local_stack, output_stack, &multiple_of_pi_to_add, &pi.max) };
    if (multiple_of_pi_to_add.numerator->sign < 0)
    {
        SWAP(pi_estimate_multiple.min, pi_estimate_multiple.max, struct Rational);
    }
    out = rational_interval_add(output_stack, local_stack, &out, &pi_estimate_multiple);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval positive_rational_float_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerDivision division =
        integer_euclidean_divide(local_stack, output_stack, a->numerator, a->denominator);
    struct Float out_get_min = { division.quotient, &zero };
    struct Integer*estimate_denominator = &one;
    while (true)
    {
        if (division.remainder->value_count == 0)
        {
            struct FloatInterval out;
            out.min = float_copy(output_stack, &out_get_min);
            out.max = out.min;
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (integer_compare(output_stack, local_stack, interval_size->denominator,
            integer_multiply(local_stack, output_stack, estimate_denominator,
                interval_size->numerator)) <= 0)
        {
            break;
        }
        division = integer_euclidean_divide(local_stack, output_stack,
            integer_double(local_stack, division.remainder), a->denominator);
        estimate_denominator = integer_double(local_stack, estimate_denominator);
        out_get_min.significand = integer_add(local_stack, division.quotient,
            integer_double(local_stack, out_get_min.significand));
        out_get_min.exponent = integer_add(local_stack, out_get_min.exponent, INT(1, -1));
    }
    struct FloatInterval out =
    { float_reduce(output_stack, local_stack, out_get_min.significand, out_get_min.exponent),
        float_reduce(output_stack, local_stack,
            integer_add(local_stack, out_get_min.significand, &one), out_get_min.exponent) };
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct FloatInterval rational_get_float_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational*a, struct Rational*interval_size)
{
    if (a->numerator->sign < 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational negative_a = rational_negate(local_stack, a);
        struct FloatInterval positive_estimate = positive_rational_float_estimate(local_stack,
            output_stack, &negative_a, interval_size);
        struct FloatInterval out =
            float_interval_negate(output_stack, local_stack, &positive_estimate);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    return positive_rational_float_estimate(output_stack, local_stack, a, interval_size);
}

#include "integer/rational/gaussian_rational/gaussian_rational.c"
#include "integer/rational/matrix/matrix.c"
#include "integer/rational/rational_interval/rational_interval.c"
#include "integer/rational/rational_polynomial/rational_polynomial.c"