#include "declarations.h"

struct NestedPolynomial*nested_polynomial_copy(struct Stack*output_stack, struct NestedPolynomial*a)
{
    return polynomial_copy(rational_polynomial_copy, output_stack, (struct Polynomial*)a);
}

bool nested_polynomial_equals(struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    return polynomial_equals(&number_field_element_operations.ring_operations,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct NestedPolynomial*nested_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    return polynomial_add(&rational_polynomial_operations.ring_operations, output_stack,
        local_stack, (struct Polynomial*)a, (struct Polynomial*)b);
}

struct NestedPolynomial*nested_polynomial_negative(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a)
{
    return polynomial_negative(&rational_polynomial_operations.ring_operations, output_stack,
        (struct Polynomial*)a);
}

struct NestedPolynomial*nested_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*minuend, struct NestedPolynomial*subtrahend)
{
    return polynomial_subtract(&number_field_element_operations.ring_operations, output_stack,
        local_stack, (struct Polynomial*)minuend, (struct Polynomial*)subtrahend);
}

struct NestedPolynomial*nested_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    return polynomial_multiply(&rational_polynomial_operations.ring_operations, output_stack,
        local_stack, (struct Polynomial*)a, (struct Polynomial*)b, 0);
}

struct NestedPolynomial*nested_polynomial_rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct RationalPolynomial*b)
{
    return polynomial_multiply_by_coefficient(&rational_polynomial_operations.ring_operations,
        output_stack, local_stack, (struct Polynomial*)a, b, 0);
}

void nested_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor)
{
    polynomial_euclidean_divide(&rational_polynomial_operations, output_stack, local_stack, out,
        (struct Polynomial*)dividend, (struct Polynomial*)divisor, 0);
}

struct NestedPolynomial*nested_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*dividend, struct NestedPolynomial*divisor)
{
    struct PolynomialDivision division;
    nested_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct NestedPolynomial*)division.remainder;
}

struct NestedPolynomial*nested_polynomial_rational_polynomial_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*dividend, struct RationalPolynomial*divisor)
{
    return polynomial_divide_by_coefficient(rational_polynomial_euclidean_quotient, output_stack,
        local_stack, (struct Polynomial*)dividend, divisor);
}

struct NestedPolynomial*nested_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a)
{
    return polynomial_derivative(rational_polynomial_integer_multiply, output_stack, local_stack,
        (struct Polynomial*)a);
}

//Correct up to multiplication by a nonzero rational.
struct RationalPolynomial*nested_polynomial_resultant(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    if (a->coefficient_count == 0 || b->coefficient_count == 0)
    {
        return (struct RationalPolynomial*)&polynomial_zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*g = &rational_polynomial_one;
    struct RationalPolynomial*h = &rational_polynomial_one;
    if (b->coefficient_count > a->coefficient_count)
    {
        POINTER_SWAP(a, b);
    }
    while (b->coefficient_count > 1)
    {
        size_t degree_size_t = a->coefficient_count - b->coefficient_count;
        struct Integer*degree = integer_from_size_t(local_stack, degree_size_t);
        struct NestedPolynomial*remainder =
            nested_polynomial_euclidean_remainder(local_stack, output_stack,
                nested_polynomial_rational_polynomial_multiply(local_stack, output_stack, a,
                    rational_polynomial_exponentiate(local_stack, output_stack,
                        b->coefficients[b->coefficient_count - 1],
                        integer_add(local_stack, degree, &one))), b);
        a = b;
        b = nested_polynomial_rational_polynomial_divide(local_stack, output_stack, remainder,
            rational_polynomial_multiply(local_stack, output_stack, g,
                rational_polynomial_exponentiate(local_stack, output_stack, h, degree)));
        g = a->coefficients[a->coefficient_count - 1];
        struct RationalPolynomial*term =
            rational_polynomial_exponentiate(local_stack, output_stack, g, degree);
        if (degree_size_t > 1)
        {
            h = rational_polynomial_euclidean_quotient(local_stack, output_stack, term,
                rational_polynomial_exponentiate(local_stack, output_stack, h,
                    integer_from_size_t(local_stack, degree_size_t - 1)));
        }
        else
        {
            h = rational_polynomial_multiply(local_stack, output_stack, term,
                rational_polynomial_exponentiate(local_stack, output_stack, h,
                    integer_from_size_t(local_stack, 1 - degree_size_t)));
        }
    }
    struct RationalPolynomial*term = rational_polynomial_exponentiate(local_stack, output_stack,
        b->coefficients[b->coefficient_count - 1],
        integer_from_size_t(local_stack, a->coefficient_count - 1));
    struct RationalPolynomial*out;
    if (a->coefficient_count > 2)
    {
        out = rational_polynomial_euclidean_quotient(output_stack, local_stack, term,
            rational_polynomial_exponentiate(local_stack, output_stack, h,
                integer_from_size_t(local_stack, a->coefficient_count - 2)));
    }
    else
    {
        out = rational_polynomial_multiply(output_stack, local_stack, term,
            rational_polynomial_exponentiate(local_stack, output_stack, h,
                integer_from_size_t(local_stack, 2 - a->coefficient_count)));
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}