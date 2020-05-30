#include "declarations.h"

struct RationalPolynomial*rational_polynomial_copy(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_copy(output_stack, &a->coefficients[i]);
    }
    return out;
}

bool rational_polynomial_equals(struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    if (a->coefficient_count != b->coefficient_count)
    {
        return false;
    }
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (!rational_equals(a->coefficients + i, b->coefficients + i))
        {
            return false;
        }
    }
    return true;
}

void rational_polynomial_trim_leading_zeroes(struct RationalPolynomial*a)
{
    for (size_t i = a->coefficient_count; i-- > 0;)
    {
        if (rational_equals(a->coefficients + i, &rational_zero))
        {
            a->coefficient_count -= 1;
        }
        else
        {
            return;
        }
    }
}

struct RationalPolynomial*rational_polynomial_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        SWAP(a, b, struct RationalPolynomial*);
    }
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] =
            rational_add(output_stack, local_stack, a->coefficients + i, b->coefficients + i);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_copy(output_stack, a->coefficients + i);
    }
    rational_polynomial_trim_leading_zeroes(out);
    return out;
}

struct RationalPolynomial*rational_polynomial_negate(struct Stack*restrict output_stack,
    struct RationalPolynomial*a)
{
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_negate(output_stack, a->coefficients + i);
    }
    return out;
}

struct RationalPolynomial*rational_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*minuend,
    struct RationalPolynomial*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = rational_polynomial_add(output_stack, local_stack, minuend,
        rational_polynomial_negate(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*rational_polynomial_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Rational*b)
{
    if (rational_equals(b, &rational_zero))
    {
        return polynomial_zero;
    }
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_multiply(output_stack, local_stack, a->coefficients + i, b);
    }
    return out;
}

struct RationalPolynomial*rational_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_zero;
    }
    struct RationalPolynomial*out = POLYNOMIAL_ALLOCATE(output_stack,
        a->coefficient_count + b->coefficient_count - 1, struct Rational);
    for (size_t i = 0; i < out->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            struct Rational product = rational_multiply(local_stack, output_stack,
                a->coefficients + i, b->coefficients + j);
            out->coefficients[i + j] =
                rational_add(output_stack, local_stack, out->coefficients + i + j, &product);
        }
    }
    rational_polynomial_trim_leading_zeroes(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void rational_polynomial_copy_coefficients(struct Stack*output_stack, struct RationalPolynomial*a)
{
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        a->coefficients[i] = rational_copy(output_stack, &a->coefficients[i]);
    }
}

struct RationalPolynomialDivision rational_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*dividend, struct RationalPolynomial*divisor)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        return (struct RationalPolynomialDivision) { polynomial_zero,
            rational_polynomial_copy(output_stack, dividend) };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomialDivision out = { POLYNOMIAL_ALLOCATE(output_stack,
        1 + dividend->coefficient_count - divisor->coefficient_count, struct Rational),
        POLYNOMIAL_ALLOCATE(output_stack, dividend->coefficient_count, struct Rational) };
    memcpy(out.remainder->coefficients, dividend->coefficients,
        dividend->coefficient_count * sizeof(struct Rational));
    struct Rational leading_coefficient_reciprocal = rational_get_reciprocal(local_stack,
        divisor->coefficients + divisor->coefficient_count - 1);
    while (out.remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out.remainder->coefficient_count;
        struct Rational quotient = rational_multiply(local_stack, output_stack,
            out.remainder->coefficients + out.remainder->coefficient_count,
            &leading_coefficient_reciprocal);
        out.quotient->coefficients[out.remainder->coefficient_count + 1 -
            divisor->coefficient_count] = quotient;
        for (size_t i = 1; i < divisor->coefficient_count; ++i)
        {
            struct Rational product = rational_multiply(local_stack, output_stack, &quotient,
                divisor->coefficients + divisor->coefficient_count - i - 1);
            out.remainder->coefficients[out.remainder->coefficient_count - i] =
                rational_subtract(local_stack, output_stack,
                    out.remainder->coefficients + out.remainder->coefficient_count - i,
                    &product);
        }
    }
    rational_polynomial_copy_coefficients(output_stack, out.quotient);
    rational_polynomial_trim_leading_zeroes(out.remainder);
    rational_polynomial_copy_coefficients(output_stack, out.remainder);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*rational_polynomial_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*base, struct Integer*exponent)
{
    if (!exponent->value_count)
    {
        return rational_polynomial_one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = rational_polynomial_one;
    while (true)
    {
        if (exponent->value[0] & 1)
        {
            out = rational_polynomial_multiply(local_stack, output_stack, out, base);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = rational_polynomial_copy(output_stack, out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = rational_polynomial_multiply(local_stack, output_stack, base, base);
    }
}

struct Integer*get_lcm_of_denominators(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Rational list[], size_t element_count)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = &one;
    for (size_t i = 0; i < element_count; ++i)
    {
        out = integer_get_lcm(local_stack, output_stack, out, list[i].denominator);
    }
    out = integer_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*rational_polynomial_to_integer_polynomial(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct Integer*multiple)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_multiply(output_stack, local_stack,
            a->coefficients[i].numerator, integer_euclidean_divide(local_stack, output_stack,
                multiple, a->coefficients[i].denominator).quotient);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct PolynomialExtendedGCDInfo rational_polynomial_get_extended_gcd(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*lcm = integer_get_lcm(local_stack, output_stack,
        get_lcm_of_denominators(local_stack, output_stack, a->coefficients, a->coefficient_count),
        get_lcm_of_denominators(local_stack, output_stack, b->coefficients, b->coefficient_count));
    struct IntegerPolynomial*integer_a =
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, a, lcm);
    struct IntegerPolynomial*integer_b =
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, b, lcm);
    struct PolynomialExtendedGCDInfo info =
        integer_polynomial_get_extended_gcd(local_stack, output_stack, integer_a, integer_b);
    struct Rational lcm_reciprocal = { INT(1, lcm->sign), lcm };
    lcm_reciprocal.denominator->sign = 1;
    struct PolynomialExtendedGCDInfo out =
    { (struct Polynomial*)rational_polynomial_rational_multiply(output_stack, local_stack,
        integer_polynomial_to_rational_polynomial(local_stack, (struct IntegerPolynomial*)info.gcd),
        &lcm_reciprocal), rational_polynomial_copy(output_stack, info.a_coefficient) };
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*rational_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = integer_polynomial_get_monic(output_stack, local_stack,
        integer_polynomial_get_gcd(local_stack, output_stack,
            rational_polynomial_to_integer_polynomial(local_stack, output_stack, a,
                get_lcm_of_denominators(local_stack, output_stack, a->coefficients,
                    a->coefficient_count)),
            rational_polynomial_to_integer_polynomial(local_stack, output_stack, b,
                get_lcm_of_denominators(local_stack, output_stack, b->coefficients,
                    b->coefficient_count))));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*rational_polynomial_get_derivative(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count - 1, struct Rational);
    struct Integer*multiplier = &zero;
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        multiplier = integer_add(local_stack, multiplier, &one);
        out->coefficients[i - 1] =
            rational_integer_multiply(output_stack, local_stack, a->coefficients + i, multiplier);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t rational_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = primitive_integer_polynomial_factor(local_stack, output_stack,
        (struct IntegerPolynomial**)out,
        rational_polynomial_get_primitive_part(local_stack, output_stack, a));
    for (size_t i = 0; i < factor_count; ++i)
    {
        out[i] = integer_polynomial_get_monic(output_stack, local_stack,
            (struct IntegerPolynomial*)out[i]);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

struct Rational rational_polynomial_evaluate_at_rational(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Rational*argument)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational out = rational_zero;
    struct Rational argument_power = rational_one;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        struct Rational product =
            rational_multiply(local_stack, output_stack, &argument_power, a->coefficients + i);
        out = rational_add(local_stack, output_stack, &out, &product);
        argument_power = rational_multiply(local_stack, output_stack, &argument_power, argument);
    }
    out = rational_copy(output_stack, &out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRational rational_polynomial_evaluate_at_gaussian_rational(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct GaussianRational*argument)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational out = gaussian_rational_zero;
    struct GaussianRational argument_power = gaussian_rational_one;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        struct GaussianRational product = gaussian_rational_rational_multiply(local_stack,
            output_stack, &argument_power, a->coefficients + i);
        out = gaussian_rational_add(local_stack, output_stack, &out, &product);
        argument_power =
            gaussian_rational_multiply(local_stack, output_stack, &argument_power, argument);
    }
    out = gaussian_rational_copy(output_stack, &out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRationalPolynomial*rational_polynomial_evaluate_at_gaussian_rational_polynomial(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct RationalPolynomial*a, struct GaussianRationalPolynomial*argument)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRationalPolynomial*out = polynomial_zero;
    struct GaussianRationalPolynomial*argument_power = gaussian_rational_polynomial_one;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out = gaussian_rational_polynomial_add(local_stack, output_stack, out,
            gaussian_rational_polynomial_rational_multiply(local_stack, output_stack,
                argument_power, a->coefficients + i));
        argument_power = gaussian_rational_polynomial_multiply(local_stack, output_stack,
            argument_power, argument);
    }
    out = gaussian_rational_polynomial_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t get_sign_variation(int8_t previous_sign, int8_t sign)
{
    int8_t out = sign - previous_sign;
    if (out < 0)
    {
        return -out;
    }
    return out;
}

size_t get_twice_cauchy_index(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*numerator,
    struct RationalPolynomial*denominator)
{
    void*local_stack_a_savepoint = local_stack_a->cursor;
    int8_t start_sign = rational_polynomial_evaluate_at_rational(local_stack_a, local_stack_b,
        denominator, &rational_zero).numerator->sign;
    size_t start_sign_variation = 0;
    int8_t end_sign = rational_polynomial_evaluate_at_rational(local_stack_a, local_stack_b,
        denominator, &rational_one).numerator->sign;
    size_t end_sign_variation = 0;
    while (numerator->coefficient_count)
    {
        int8_t next_start_sign = rational_polynomial_evaluate_at_rational(local_stack_a,
            local_stack_b, numerator, &rational_zero).numerator->sign;
        start_sign_variation += get_sign_variation(start_sign, next_start_sign);
        start_sign = next_start_sign;
        int8_t next_end_sign = rational_polynomial_evaluate_at_rational(local_stack_a,
            local_stack_b, numerator, &rational_one).numerator->sign;
        end_sign_variation += get_sign_variation(end_sign, next_end_sign);
        end_sign = next_end_sign;
        struct RationalPolynomial*remainder = rational_polynomial_negate(local_stack_a,
            rational_polynomial_euclidean_divide(local_stack_a, local_stack_b, denominator,
                numerator).remainder);
        denominator = numerator;
        numerator = remainder;
    }
    local_stack_a->cursor = local_stack_a_savepoint;
    return start_sign_variation - end_sign_variation;
}

void rational_polynomial_parameterize_over_segment(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial**out_real,
    struct RationalPolynomial**out_imaginary, struct RationalPolynomial*a,
    struct GaussianRational*endpoint_a, struct GaussianRational*endpoint_b)
{
    struct GaussianRational difference =
        gaussian_rational_subtract(local_stack, output_stack, endpoint_b, endpoint_a);
    POLY(argument, struct GaussianRationalPolynomial, struct GaussianRational, 2, *endpoint_a,
        difference);
    struct GaussianRationalPolynomial*a_at_edge =
        rational_polynomial_evaluate_at_gaussian_rational_polynomial(local_stack, output_stack, a,
            &argument.p);
    *out_real = POLYNOMIAL_ALLOCATE(output_stack, a_at_edge->coefficient_count, struct Rational);
    *out_imaginary =
        POLYNOMIAL_ALLOCATE(output_stack, a_at_edge->coefficient_count, struct Rational);
    for (size_t i = 0; i < a_at_edge->coefficient_count; ++i)
    {
        (*out_real)->coefficients[i] =
            rational_copy(output_stack, &a_at_edge->coefficients[i].real);
        (*out_imaginary)->coefficients[i] =
            rational_copy(output_stack, &a_at_edge->coefficients[i].imaginary);
    }
    rational_polynomial_trim_leading_zeroes(*out_real);
    rational_polynomial_trim_leading_zeroes(*out_imaginary);
}

size_t rational_polynomial_count_roots_over_segment(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*a,
    struct GaussianRational*endpoint_a, struct GaussianRational*endpoint_b)
{
    void*local_stack_a_savepoint = local_stack_a->cursor;
    struct RationalPolynomial*parameterization_real_part;
    struct RationalPolynomial*parameterization_imaginary_part;
    rational_polynomial_parameterize_over_segment(local_stack_a, local_stack_b,
        &parameterization_real_part, &parameterization_imaginary_part, a, endpoint_a, endpoint_b);
    struct RationalPolynomial*sturm_chain_first_element = rational_polynomial_get_gcd(local_stack_a,
        local_stack_b, parameterization_real_part, parameterization_imaginary_part);
    size_t out = get_twice_cauchy_index(local_stack_a, local_stack_b,
        rational_polynomial_get_derivative(local_stack_a, local_stack_b, sturm_chain_first_element),
        sturm_chain_first_element);
    struct GaussianRational evaluation =
        rational_polynomial_evaluate_at_gaussian_rational(local_stack_a, local_stack_b, a,
            endpoint_a);
    if (gaussian_rational_equals(&gaussian_rational_zero, &evaluation))
    {
        ++out;
    }
    evaluation = rational_polynomial_evaluate_at_gaussian_rational(local_stack_a, local_stack_b, a,
        endpoint_b);
    if (gaussian_rational_equals(&gaussian_rational_zero, &evaluation))
    {
        ++out;
    }
    local_stack_a->cursor = local_stack_a_savepoint;
    ASSERT(out % 2 == 0,
        "rational_polynomial_count_roots_over_segment found non-integer number of roots.");
    return out / 2;
}

//Inconclusive if the rectangle is nondegerate and has a zero on a vertex; returns
//a->coefficient_count in that case.
size_t rational_polynomial_count_roots_in_rectangle(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*a, struct RationalInterval*real,
    struct RationalInterval*imaginary)
{
    struct GaussianRational corners[5] = { { real->min, imaginary->min },
        { real->max, imaginary->min }, { real->max, imaginary->max }, { real->min, imaginary->max },
        { real->min, imaginary->min } };
    if (rational_equals(&imaginary->min, &imaginary->max))
    {
        if (rational_equals(&real->min, &real->max))
        {
            struct GaussianRational evaluation =
                rational_polynomial_evaluate_at_gaussian_rational(local_stack_a, local_stack_b, a,
                    corners);
            if (gaussian_rational_equals(&gaussian_rational_zero, &evaluation))
            {
                return 1;
            }
            return 0;
        }
        return rational_polynomial_count_roots_over_segment(local_stack_a, local_stack_b, a,
            corners, corners + 1);
    }
    if (rational_equals(&real->min, &real->max))
    {
        return rational_polynomial_count_roots_over_segment(local_stack_a, local_stack_b, a,
            corners + 1, corners + 2);
    }
    for (size_t i = 0; i < 4; ++i)
    {
        struct GaussianRational evaluation =
            rational_polynomial_evaluate_at_gaussian_rational(local_stack_a, local_stack_b, a,
                corners + i);
        if (gaussian_rational_equals(&gaussian_rational_zero, &evaluation))
        {
            return a->coefficient_count;
        }
    }
    size_t out = 0;
    for (size_t i = 0; i < 4; ++i)
    {
        void*local_stack_a_savepoint = local_stack_a->cursor;
        struct RationalPolynomial*parameterization_real_part;
        struct RationalPolynomial*parameterization_imaginary_part;
        rational_polynomial_parameterize_over_segment(local_stack_a, local_stack_b,
            &parameterization_real_part, &parameterization_imaginary_part, a, corners + i,
            corners + i + 1);
        out += get_twice_cauchy_index(local_stack_a, local_stack_b, parameterization_real_part,
            parameterization_imaginary_part);
        struct RationalPolynomial*sturm_chain_first_element =
            rational_polynomial_get_gcd(local_stack_a, local_stack_b, parameterization_real_part,
                parameterization_imaginary_part);
        out += get_twice_cauchy_index(local_stack_a, local_stack_b,
            rational_polynomial_get_derivative(local_stack_a, local_stack_b,
                sturm_chain_first_element),
            sturm_chain_first_element);
        local_stack_a->cursor = local_stack_a_savepoint;
    }
    ASSERT(out % 4 == 0,
        "rational_polynomial_count_roots_in_rectangle found non-integer number of roots.");
    return out / 4;
}

size_t rational_polynomial_count_roots_in_region(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct RationalPolynomial*a, struct Region*region)
{
    void*local_stack_a_savepoint = local_stack_a->cursor;
    struct RationalInterval real_interval =
        float_interval_to_rational_interval(local_stack_a, local_stack_b, &region->real_interval);
    struct RationalInterval imaginary_interval = float_interval_to_rational_interval(local_stack_a,
        local_stack_b, &region->imaginary_interval);
    size_t out = rational_polynomial_count_roots_in_rectangle(local_stack_a, local_stack_b, a,
        &real_interval, &imaginary_interval);
    local_stack_a->cursor = local_stack_a_savepoint;
    return out;
}

//In this case, interval_size is only passed through to internal calls to rational_float_estimate
//rather than indicating that the intervals in out will be no larger than interval_size. It's still
//the case that the sizes of the intervals in out approach 0 as interval_size and the sizes of the
//intervals in argument approach 0. This is all that turns out to be necessary in this program, so
//no more direct means of controlling the accuracy of out are provided.
struct Region rational_polynomial_evaluate_at_region_estimate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct Region*argument,
    struct Rational*interval_size)
{
    if (a->coefficient_count == 0)
    {
        return (struct Region) { (struct FloatInterval) { float_zero, float_zero },
            (struct FloatInterval) { float_zero, float_zero } };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval*coefficients =
        ARRAY_ALLOCATE(local_stack, a->coefficient_count, struct FloatInterval);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        coefficients[i] = rational_get_float_estimate(local_stack, output_stack,
            a->coefficients + i, interval_size);
    }
    struct Region out = { coefficients[0].min, coefficients[0].max, float_zero, float_zero };
    struct Region argument_power = { { float_one, float_one }, { float_zero, float_zero } };
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        struct FloatInterval product_a = float_interval_multiply(local_stack, output_stack,
            &argument_power.real_interval, &argument->real_interval);
        struct FloatInterval product_b = float_interval_multiply(local_stack, output_stack,
            &argument_power.imaginary_interval, &argument->imaginary_interval);
        struct FloatInterval new_argument_power_real_part =
            float_interval_subtract(local_stack, output_stack, &product_a, &product_b);
        product_a = float_interval_multiply(local_stack, output_stack,
            &argument_power.real_interval, &argument->imaginary_interval);
        product_b = float_interval_multiply(local_stack, output_stack,
            &argument_power.imaginary_interval, &argument->real_interval);
        argument_power.imaginary_interval =
            float_interval_add(local_stack, output_stack, &product_a, &product_b);
        argument_power.real_interval = new_argument_power_real_part;
        product_a = float_interval_multiply(local_stack, output_stack,
            &argument_power.real_interval, coefficients + i);
        out.real_interval =
            float_interval_add(local_stack, output_stack, &product_a, &out.real_interval);
        product_a = float_interval_multiply(local_stack, output_stack,
            &argument_power.imaginary_interval, coefficients + i);
        out.imaginary_interval =
            float_interval_add(local_stack, output_stack, &product_a, &out.imaginary_interval);
    }
    out.real_interval = float_interval_copy(output_stack, &out.real_interval);
    out.imaginary_interval = float_interval_copy(output_stack, &out.imaginary_interval);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*rational_polynomial_get_primitive_part(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_get_primitive_part(output_stack, local_stack,
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, a, 
            get_lcm_of_denominators(local_stack, output_stack, a->coefficients,
                a->coefficient_count)));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (a->coefficients[i].numerator->value_count)
        {
            out->coefficients[i] = POLYNOMIAL_ALLOCATE(output_stack, 1, struct Rational);
            out->coefficients[i]->coefficients[0] =
                rational_copy(output_stack, a->coefficients + i);
        }
        else
        {
            out->coefficients[i] = polynomial_zero;
        }
    }
    return out;
}

#include "integer/rational/rational_polynomial/nested_polynomial/nested_polynomial.c"
#include "integer/rational/rational_polynomial/number_field_element/number_field_element.c"