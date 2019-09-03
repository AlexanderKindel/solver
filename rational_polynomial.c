#include "declarations.h"

struct RationalPolynomial*rational_polynomial_copy(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    return polynomial_copy(rational_copy, output_stack, (struct Polynomial*)a);
}

bool rational_polynomial_equals(struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_equals(&rational_operations.ring_operations, (struct Polynomial*)a,
        (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_add(&rational_operations.ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_negative(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    return polynomial_negative(&rational_operations.ring_operations, output_stack,
        (struct Polynomial*)a);
}

struct RationalPolynomial*rational_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_subtract(&rational_operations.ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct RationalPolynomial*rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_multiply(&rational_operations.ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b, 0);
}

struct RationalPolynomial*rational_polynomial_integer_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Integer*b)
{
    return rational_polynomial_rational_multiply(output_stack, local_stack, a,
        &(struct Rational){ b, &one });
}

struct RationalPolynomial*rational_polynomial_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Rational*b)
{
    return polynomial_multiply_by_coefficient(&rational_operations.ring_operations, output_stack,
        local_stack, (struct Polynomial*)a, b, 0);
}

void rational_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor)
{
    field_polynomial_euclidean_divide(&rational_operations, output_stack, local_stack, out,
        (struct Polynomial*)dividend, (struct Polynomial*)divisor, 0);
}

struct RationalPolynomial*rational_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend, struct RationalPolynomial*divisor)
{
    struct PolynomialDivision division;
    rational_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct RationalPolynomial*)division.quotient;
}

struct RationalPolynomial*rational_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend, struct RationalPolynomial*divisor)
{
    struct PolynomialDivision division;
    rational_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct RationalPolynomial*)division.remainder;
}

struct Integer*denominator_lcm(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*list[], size_t element_count)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = &one;
    for (size_t i = 0; i < element_count; ++i)
    {
        out = integer_lcm(local_stack, output_stack, out, list[i]->denominator);
    }
    out = integer_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*rational_polynomial_to_integer_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Integer*multiple)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_multiply(output_stack, local_stack,
            a->coefficients[i]->numerator, integer_euclidean_quotient(local_stack, output_stack,
                multiple, a->coefficients[i]->denominator));
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*rational_polynomial_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*base, struct Integer*exponent)
{
    return generic_exponentiate(&rational_polynomial_operations.ring_operations, output_stack,
        local_stack, base, exponent, 0);
}

void rational_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialExtendedGCDInfo*out, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*lcm = integer_lcm(local_stack, output_stack,
        denominator_lcm(local_stack, output_stack, a->coefficients, a->coefficient_count),
        denominator_lcm(local_stack, output_stack, b->coefficients, b->coefficient_count));
    struct IntegerPolynomial*integer_a =
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, a, lcm);
    struct IntegerPolynomial*integer_b =
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, b, lcm);
    struct PolynomialExtendedGCDInfo info;
    integer_polynomial_extended_gcd(local_stack, output_stack, &info, integer_a, integer_b);
    struct Rational lcm_reciprocal = { &(struct Integer) { 1, lcm->sign, 1 }, lcm };
    lcm_reciprocal.denominator->sign = 1;
    out->gcd = (struct Polynomial*)rational_polynomial_rational_multiply(output_stack, local_stack,
        integer_polynomial_to_rational_polynomial(local_stack, (struct IntegerPolynomial*)info.gcd),
        &lcm_reciprocal);
    out->a_coefficient = rational_polynomial_copy(output_stack, info.a_coefficient);
    local_stack->cursor = local_stack_savepoint;
}

struct RationalPolynomial*rational_polynomial_gcd(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct PolynomialExtendedGCDInfo info;
    integer_polynomial_extended_gcd(local_stack, output_stack, &info,
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, a,
            denominator_lcm(local_stack, output_stack, a->coefficients, a->coefficient_count)),
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, b,
            denominator_lcm(local_stack, output_stack, b->coefficients, b->coefficient_count)));
    struct RationalPolynomial*out =
        integer_polynomial_to_monic(output_stack, local_stack, (struct IntegerPolynomial*)info.gcd);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t rational_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*a, struct RationalPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = primitive_integer_polynomial_factor(local_stack, output_stack,
        rational_polynomial_primitive_part(local_stack, output_stack, a),
        (struct IntegerPolynomial**)out);
    for (size_t i = 0; i < factor_count; ++i)
    {
        out[i] = integer_polynomial_to_monic(output_stack, local_stack,
            (struct IntegerPolynomial*)out[i]);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

struct RationalPolynomial*rational_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a)
{
    return polynomial_derivative(rational_integer_multiply, output_stack, local_stack,
        (struct Polynomial*)a);
}

void*rational_polynomial_evaluate(struct RingOperations*argument_operations,
    void*(argument_rational_multiply)(struct Stack*, struct Stack*, void*, struct Rational*),
    struct Stack*output_stack, struct Stack*local_stack, struct RationalPolynomial*a, void*argument)
{
    void*local_stack_savepoint = local_stack->cursor;
    void*out = argument_operations->additive_identity;
    void*argument_power = argument_operations->multiplicative_identity;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out = argument_operations->add(local_stack, output_stack, out,
            argument_rational_multiply(local_stack, output_stack, argument_power,
                a->coefficients[i]), 0);
        argument_power =
            argument_operations->multiply(local_stack, output_stack, argument_power, argument, 0);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational*rational_polynomial_evaluate_at_rational(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Rational*argument)
{
    return rational_polynomial_evaluate(&rational_operations.ring_operations, rational_multiply,
        output_stack, local_stack, a, argument);
}

struct GaussianRational*rational_polynomial_evaluate_at_gaussian_rational(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct GaussianRational*argument)
{
    return rational_polynomial_evaluate(&gaussian_rational_operations,
        gaussian_rational_rational_multiply, output_stack, local_stack, a, argument);
}

size_t sign_variation(int8_t previous_sign, int8_t sign)
{
    int8_t out = sign - previous_sign;
    if (out < 0)
    {
        return -out;
    }
    return out;
}

size_t twice_cauchy_index(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalPolynomial*numerator, struct RationalPolynomial*denominator)
{
    void*stack_a_savepoint = stack_a->cursor;
    int8_t start_sign = rational_polynomial_evaluate_at_rational(stack_a, stack_b, denominator,
        &rational_zero)->numerator->sign;
    size_t start_sign_variation = 0;
    int8_t end_sign = rational_polynomial_evaluate_at_rational(stack_a, stack_b, denominator,
        &rational_one)->numerator->sign;
    size_t end_sign_variation = 0;
    while (numerator->coefficient_count)
    {
        int8_t next_start_sign = rational_polynomial_evaluate_at_rational(stack_a, stack_b,
            numerator, &rational_zero)->numerator->sign;
        start_sign_variation += sign_variation(start_sign, next_start_sign);
        start_sign = next_start_sign;
        int8_t next_end_sign = rational_polynomial_evaluate_at_rational(stack_a, stack_b,
            numerator, &rational_one)->numerator->sign;
        end_sign_variation += sign_variation(end_sign, next_end_sign);
        end_sign = next_end_sign;
        struct RationalPolynomial*remainder = rational_polynomial_negative(stack_a, 
            rational_polynomial_euclidean_remainder(stack_a, stack_b, denominator, numerator));
        denominator = numerator;
        numerator = remainder;
    }
    stack_a->cursor = stack_a_savepoint;
    return start_sign_variation - end_sign_variation;
}

void rational_polynomial_parameterize_over_segment(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial**out_real,
    struct RationalPolynomial**out_imaginary, struct RationalPolynomial*a,
    struct GaussianRational*endpoint_a, struct GaussianRational*endpoint_b)
{
    struct GaussianRationalPolynomial*a_at_edge =
        rational_polynomial_evaluate(&gaussian_rational_polynomial_operations,
            gaussian_rational_polynomial_rational_multiply, local_stack, output_stack, a,
            &(struct GaussianRationalPolynomial){2, endpoint_a,
                gaussian_rational_subtract(local_stack, output_stack, endpoint_b, endpoint_a)});
    *out_real = polynomial_allocate(output_stack, a_at_edge->coefficient_count);
    *out_imaginary = polynomial_allocate(output_stack, a_at_edge->coefficient_count);
    for (size_t i = 0; i < a_at_edge->coefficient_count; ++i)
    {
        (*out_real)->coefficients[i] =
            rational_copy(output_stack, a_at_edge->coefficients[i]->real);
        (*out_imaginary)->coefficients[i] =
            rational_copy(output_stack, a_at_edge->coefficients[i]->imaginary);
    }
    polynomial_trim_leading_zeroes(&rational_operations.ring_operations,
        (struct Polynomial*)(*out_real));
    polynomial_trim_leading_zeroes(&rational_operations.ring_operations,
        (struct Polynomial*)(*out_imaginary));
}

size_t rational_polynomial_root_count_over_segment(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalPolynomial*a, struct GaussianRational*endpoint_a,
    struct GaussianRational*endpoint_b)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct RationalPolynomial*parameterization_real_part;
    struct RationalPolynomial*parameterization_imaginary_part;
    rational_polynomial_parameterize_over_segment(stack_a, stack_b, &parameterization_real_part,
        &parameterization_imaginary_part, a, endpoint_a, endpoint_b);
    struct RationalPolynomial*sturm_chain_first_element = rational_polynomial_gcd(stack_a, stack_b,
        parameterization_real_part, parameterization_imaginary_part);
    size_t out = twice_cauchy_index(stack_a, stack_b,
        rational_polynomial_derivative(stack_a, stack_b, sturm_chain_first_element),
        sturm_chain_first_element);
    if (gaussian_rational_equals(&gaussian_rational_zero,
        rational_polynomial_evaluate_at_gaussian_rational(stack_a, stack_b, a, endpoint_a)))
    {
        ++out;
    }
    if (gaussian_rational_equals(&gaussian_rational_zero,
        rational_polynomial_evaluate_at_gaussian_rational(stack_a, stack_b, a, endpoint_b)))
    {
        ++out;
    }
    stack_a->cursor = stack_a_savepoint;
    ASSERT(out % 2 == 0,
        "rational_polynomial_root_count_over_segment found fractional number of roots.");
    return out / 2;
}

//Inconclusive if the rectangle is nondegerate and has a zero on a vertex; returns
//a->coefficient_count in that case.
size_t rational_polynomial_root_count_in_rectangle(struct Stack*stack_a, struct Stack*stack_b,
    struct RationalPolynomial*a, struct RationalInterval*real, struct RationalInterval*imaginary)
{
    struct GaussianRational corners[5] = { { real->min, imaginary->min },
        { real->max, imaginary->min }, { real->max, imaginary->max }, { real->min, imaginary->max },
        { real->min, imaginary->min } };
    if (rational_equals(imaginary->min, imaginary->max))
    {
        if (rational_equals(real->min, real->max))
        {
            if (gaussian_rational_equals(&gaussian_rational_zero,
                rational_polynomial_evaluate_at_gaussian_rational(stack_a, stack_b, a, corners)))
            {
                return 1;
            }
            return 0;
        }
        return rational_polynomial_root_count_over_segment(stack_a, stack_b, a, corners,
            corners + 1);
    }
    if (rational_equals(real->min, real->max))
    {
        return rational_polynomial_root_count_over_segment(stack_a, stack_b, a, corners + 1,
            corners + 2);
    }
    for (size_t i = 0; i < 4; ++i)
    {
        if (gaussian_rational_equals(&gaussian_rational_zero,
            rational_polynomial_evaluate_at_gaussian_rational(stack_a, stack_b, a, corners + i)))
        {
            return a->coefficient_count;
        }
    }
    size_t out = 0;
    for (size_t i = 0; i < 4; ++i)
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct RationalPolynomial*parameterization_real_part;
        struct RationalPolynomial*parameterization_imaginary_part;
        rational_polynomial_parameterize_over_segment(stack_a, stack_b, &parameterization_real_part,
            &parameterization_imaginary_part, a, corners + i, corners + i + 1);
        out += twice_cauchy_index(stack_a, stack_b, parameterization_real_part,
            parameterization_imaginary_part);
        struct RationalPolynomial*sturm_chain_first_element = rational_polynomial_gcd(stack_a,
            stack_b, parameterization_real_part, parameterization_imaginary_part);
        out += twice_cauchy_index(stack_a, stack_b,
            rational_polynomial_derivative(stack_a, stack_b, sturm_chain_first_element),
            sturm_chain_first_element);
        stack_a->cursor = stack_a_savepoint;
    }
    ASSERT(out % 4 == 0,
        "rational_polynomial_root_count_in_rectangle found fractional number of roots.");
    return out / 4;
}

//In this case, interval_size is only passed through to internal calls to rational_float_estimate
//rather than indicating that the intervals in out will be no larger than interval_size. It's still
//the case that the sizes of the intervals in out converge to 0 as interval_size and the sizes of
//the intervals in argument converge to 0. This is all that turns out to be necessary in this
//program, so no more direct means of controlling the accuracy of out are provided.
void rational_polynomial_evaluate_at_rectangular_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct RationalPolynomial*a,
    struct RectangularEstimate*argument, struct Rational*interval_size)
{
    out->imaginary_part_estimate = ALLOCATE(output_stack, struct FloatInterval);
    if (a->coefficient_count == 0)
    {
        out->real_part_estimate = ALLOCATE(output_stack, struct FloatInterval);
        out->real_part_estimate->min = &float_zero;
        out->real_part_estimate->max = &float_zero;
        out->imaginary_part_estimate->min = &float_zero;
        out->imaginary_part_estimate->max = &float_zero;
        return;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval**coefficients =
        ARRAY_ALLOCATE(local_stack, a->coefficient_count, struct FloatInterval*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        coefficients[i] =
            rational_float_estimate(local_stack, output_stack, a->coefficients[i], interval_size);
    }
    out->real_part_estimate = ALLOCATE(output_stack, struct FloatInterval);
    out->real_part_estimate->min = coefficients[0]->min;
    out->real_part_estimate->max = coefficients[0]->max;
    out->imaginary_part_estimate->min = &float_zero;
    out->imaginary_part_estimate->max = &float_zero;
    struct RectangularEstimate argument_power = { &(struct FloatInterval){&float_one, &float_one},
        &(struct FloatInterval){&float_zero, &float_zero} };
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        struct FloatInterval*new_argument_power_real_part =
            float_interval_subtract(local_stack, output_stack,
                float_interval_multiply(local_stack, output_stack,
                    argument_power.real_part_estimate, argument->real_part_estimate),
                float_interval_multiply(local_stack, output_stack,
                    argument_power.imaginary_part_estimate, argument->imaginary_part_estimate));
        struct FloatInterval*new_argument_power_imaginary_part =
            float_interval_add(local_stack, output_stack,
                float_interval_multiply(local_stack, output_stack,
                    argument_power.real_part_estimate, argument->imaginary_part_estimate),
                float_interval_multiply(local_stack, output_stack,
                    argument_power.imaginary_part_estimate, argument->real_part_estimate));
        argument_power.real_part_estimate = new_argument_power_real_part;
        argument_power.imaginary_part_estimate = new_argument_power_imaginary_part;
        out->real_part_estimate = float_interval_add(local_stack, output_stack,
            float_interval_multiply(local_stack, output_stack, argument_power.real_part_estimate,
                coefficients[i]), out->real_part_estimate);
        out->imaginary_part_estimate = float_interval_add(local_stack, output_stack,
            float_interval_multiply(local_stack, output_stack,
                argument_power.imaginary_part_estimate, coefficients[i]),
            out->imaginary_part_estimate);
    }
    out->real_part_estimate = float_interval_copy(output_stack, out->real_part_estimate);
    out->imaginary_part_estimate = float_interval_copy(output_stack, out->imaginary_part_estimate);
    local_stack->cursor = local_stack_savepoint;
}

struct IntegerPolynomial*rational_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_primitive_part(output_stack, local_stack,
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, a, 
            denominator_lcm(local_stack, output_stack, a->coefficients, a->coefficient_count)));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    struct NestedPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (a->coefficients[i]->numerator->value_count)
        {
            out->coefficients[i] = polynomial_allocate(output_stack, 1);
            out->coefficients[i]->coefficients[0] = rational_copy(output_stack, a->coefficients[i]);
        }
        else
        {
            out->coefficients[i] = (struct RationalPolynomial*)&polynomial_zero;
        }
    }
    return out;
}