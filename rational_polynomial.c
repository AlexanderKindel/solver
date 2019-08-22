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

struct Integer*sturm_chain_sign_variation(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*first_element, struct RationalPolynomial*second_element)
{
    void*local_stack_savepoint = local_stack->cursor;
    int8_t start_sign = rational_polynomial_evaluate_at_rational(local_stack, output_stack,
        first_element, &rational_zero)->numerator->sign;
    size_t start_sign_change_count = 0;
    int8_t end_sign = rational_polynomial_evaluate_at_rational(local_stack, output_stack,
        first_element, &rational_one)->numerator->sign;
    size_t end_sign_change_count = 0;
    while (second_element->coefficient_count > 1)
    {
        struct Rational*second_element_at_edge_start =
            rational_polynomial_evaluate_at_rational(local_stack, output_stack, second_element,
                &rational_zero);
        int8_t sign_change = second_element_at_edge_start->numerator->sign - start_sign;
        start_sign_change_count += sign_change * sign_change;
        start_sign = second_element_at_edge_start->numerator->sign;
        struct Rational*second_element_at_edge_end =
            rational_polynomial_evaluate_at_rational(local_stack, output_stack, second_element,
                &rational_one);
        sign_change = second_element_at_edge_end->numerator->sign - end_sign;
        end_sign_change_count += sign_change * sign_change;
        end_sign = second_element_at_edge_end->numerator->sign;
        struct RationalPolynomial*remainder = rational_polynomial_negative(local_stack, 
            rational_polynomial_euclidean_remainder(local_stack, output_stack, first_element,
                second_element));
        first_element = second_element;
        second_element = remainder;
    }
    struct Integer*out = integer_subtract(output_stack, local_stack,
        integer_from_size_t(local_stack, start_sign_change_count),
        integer_from_size_t(local_stack, end_sign_change_count));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Returns a->coefficient_count if a has a zero on the rectangle boundary.
size_t rational_polynomial_root_count_in_rectangle(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalInterval*real,
    struct RationalInterval*imaginary)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational corners[5] = { { real->min, imaginary->min },
        { real->max, imaginary->min }, { real->max, imaginary->max }, { real->min, imaginary->max },
        { real->min, imaginary->min } };
    struct GaussianRationalPolynomial*edges[4];
    struct Integer*sign_variation_total = &zero;
    for (size_t i = 0; i < 4; ++i)
    {
        edges[i] = polynomial_allocate(local_stack, 2);
        edges[i]->coefficients[0] = &corners[i];
        edges[i]->coefficients[1] =
            gaussian_rational_subtract(local_stack, output_stack, &corners[i + 1], &corners[i]);
        struct GaussianRationalPolynomial*a_at_edge =
            rational_polynomial_evaluate(&gaussian_rational_polynomial_operations,
                gaussian_rational_polynomial_rational_multiply, local_stack, output_stack, a,
                edges[i]);
        struct RationalPolynomial*real_part =
            polynomial_allocate(local_stack, a->coefficient_count);
        struct RationalPolynomial*imaginary_part =
            polynomial_allocate(local_stack, a->coefficient_count);
        for (size_t i = 0; i < a->coefficient_count; ++i)
        {
            real_part->coefficients[i] =
                rational_copy(local_stack, a_at_edge->coefficients[i]->real);
            imaginary_part->coefficients[i] =
                rational_copy(local_stack, a_at_edge->coefficients[i]->imaginary);
        }
        polynomial_trim_leading_zeroes(&rational_operations.ring_operations,
            (struct Polynomial*)real_part);
        polynomial_trim_leading_zeroes(&rational_operations.ring_operations,
            (struct Polynomial*)imaginary_part);
        struct RationalPolynomial*sturm_chain_first_element =
            rational_polynomial_gcd(local_stack, output_stack, real_part, imaginary_part);
        if (sturm_chain_sign_variation(local_stack, output_stack, sturm_chain_first_element,
            rational_polynomial_derivative(local_stack, output_stack,
                sturm_chain_first_element))->sign)
        {
            local_stack->cursor = local_stack_savepoint;
            return a->coefficient_count;
        }
        sign_variation_total = integer_add(local_stack, sign_variation_total,
            sturm_chain_sign_variation(local_stack, output_stack, real_part, imaginary_part));
    }
    sign_variation_total->sign = -sign_variation_total->sign;
    size_t out = integer_to_size_t(integer_half(output_stack, sign_variation_total));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void rational_polynomial_evaluate_at_rectangular_estimate(struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct RationalPolynomial*a,
    struct RectangularEstimate*argument)
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
    struct Rational*coefficient_interval_size = rational_integer_divide(local_stack, output_stack,
        float_to_rational(local_stack, output_stack,
            float_subtract(local_stack, output_stack, argument->real_part_estimate->max,
                argument->real_part_estimate->min)),
        integer_from_size_t(local_stack, a->coefficient_count));
    struct FloatInterval**coefficients =
        ARRAY_ALLOCATE(local_stack, a->coefficient_count, struct FloatInterval*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        coefficients[i] = rational_float_estimate(local_stack, output_stack, a->coefficients[i],
            coefficient_interval_size);
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