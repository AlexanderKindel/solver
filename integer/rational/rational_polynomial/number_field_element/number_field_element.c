#include "declarations.h"

struct RationalPolynomial*number_field_element_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = rational_polynomial_euclidean_divide(output_stack, local_stack,
        rational_polynomial_multiply(local_stack, output_stack, a, b),
        generator_minimal_polynomial).remainder;
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct PolynomialExtendedGCDInfo info = rational_polynomial_get_extended_gcd(local_stack,
        output_stack, a, generator_minimal_polynomial);
    struct Rational reciprocal =
        rational_get_reciprocal(local_stack, ((struct RationalPolynomial*)info.gcd)->coefficients);
    struct RationalPolynomial*out = rational_polynomial_rational_multiply(output_stack, local_stack,
        (struct RationalPolynomial*)info.a_coefficient, &reciprocal);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = number_field_element_multiply(output_stack, local_stack,
        dividend, number_field_element_get_reciprocal(local_stack, output_stack, divisor,
            generator_minimal_polynomial), generator_minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*base, struct Integer*exponent,
    struct RationalPolynomial*generator_minimal_polynomial)
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
            out = number_field_element_multiply(local_stack, output_stack, out, base,
                generator_minimal_polynomial);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = rational_polynomial_copy(output_stack, out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = number_field_element_multiply(local_stack, output_stack, base, base,
            generator_minimal_polynomial);
    }
}

#include "integer/rational/rational_polynomial/number_field_element/number_field_polynomial/number_field_polynomial.c"