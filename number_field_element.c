#include "declarations.h"

struct RationalPolynomial*number_field_element_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out =
        rational_polynomial_euclidean_remainder(output_stack, local_stack,
            rational_polynomial_multiply(local_stack, output_stack, a, b),
            generator_minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_reciprocal(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct PolynomialExtendedGCDInfo info;
    rational_polynomial_extended_gcd(local_stack, output_stack, &info, a,
        generator_minimal_polynomial);
    struct RationalPolynomial*out =
        rational_polynomial_rational_multiply(output_stack, local_stack,
        (struct RationalPolynomial*)info.a_coefficient,
            rational_reciprocal(local_stack, output_stack, info.gcd->coefficients[0]));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*dividend, struct RationalPolynomial*divisor,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = number_field_element_multiply(output_stack, local_stack,
        dividend, number_field_element_reciprocal(local_stack, output_stack, divisor,
            generator_minimal_polynomial), generator_minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*base, struct Integer*exponent,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    return generic_exponentiate(&number_field_element_operations.ring_operations, output_stack,
        local_stack, base, exponent, generator_minimal_polynomial);
}