#include "declarations.h"

struct RationalPolynomial*rational_polynomial_copy_to_stack(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    return polynomial_copy(rational_copy_to_stack, output_stack, (struct Polynomial*)a);
}

struct RationalPolynomial*rational_polynomial_copy_to_pool(struct PoolSet*pool_set,
    struct RationalPolynomial*a)
{
    struct RationalPolynomial*out = pool_polynomial_allocate(pool_set, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_copy_to_pool(pool_set, a->coefficients[i]);
    }
    return out;
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
    return polynomial_multiply_by_coefficient(&integer_operations.ring_operations, output_stack,
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
        out = integer_multiply(local_stack, output_stack, out, list[i]->denominator);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*rational_polynomial_to_integer_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct Integer*multiple)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = stack_polynomial_allocate(local_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_multiply(local_stack, output_stack,
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
    struct Integer*content = integer_gcd(local_stack, output_stack,
        integer_polynomial_content(local_stack, output_stack, integer_a),
        integer_polynomial_content(local_stack, output_stack, integer_b));
    struct PolynomialExtendedGCDInfo info;
    integer_polynomial_extended_gcd(local_stack, output_stack, &info,
        integer_polynomial_integer_divide(local_stack, output_stack, integer_a, content),
        integer_polynomial_integer_divide(local_stack, output_stack, integer_b, content));
    out->a_coefficient =
        stack_polynomial_allocate(output_stack, info.a_coefficient->coefficient_count);
    for (size_t i = 0; i < info.a_coefficient->coefficient_count; ++i)
    {
        out->a_coefficient->coefficients[i] = STACK_SLOT_ALLOCATE(output_stack, struct Rational);
        ((struct Rational*)out->a_coefficient->coefficients[i])->numerator =
            integer_copy_to_stack(output_stack, info.a_coefficient->coefficients[i]);
        ((struct Rational*)out->a_coefficient->coefficients[i])->denominator =
            stack_integer_initialize(output_stack, 1, 1);
    }
    out->gcd = stack_polynomial_allocate(output_stack, info.gcd->coefficient_count);
    for (size_t i = 0; i < info.gcd->coefficient_count; ++i)
    {
        out->gcd->coefficients[i] = rational_reduced(output_stack, local_stack,
            integer_multiply(local_stack, output_stack, content, info.gcd->coefficients[i]), lcm);
    }
    local_stack->cursor = local_stack_savepoint;
}

struct RationalPolynomial*rational_polynomial_gcd(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b)
{
    struct PolynomialExtendedGCDInfo info;
    rational_polynomial_extended_gcd(output_stack, local_stack, &info, a, b);
    return (struct RationalPolynomial*)info.gcd;
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
        out[i] = integer_polynomial_to_monic(local_stack, output_stack,
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

struct IntegerPolynomial*rational_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct RationalPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*lcm =
        denominator_lcm(local_stack, output_stack, a->coefficients, a->coefficient_count);
    struct IntegerPolynomial*out = integer_polynomial_primitive_part(output_stack, local_stack,
        rational_polynomial_to_integer_polynomial(local_stack, output_stack, a, lcm));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*rational_polynomial_to_nested_polynomial(struct Stack*output_stack,
    struct RationalPolynomial*a)
{
    struct NestedPolynomial*out = stack_polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = stack_polynomial_allocate(output_stack, 0);
        out->coefficients[i]->coefficients[0] =
            rational_copy_to_stack(output_stack, a->coefficients[i]);
    }
    return out;
}