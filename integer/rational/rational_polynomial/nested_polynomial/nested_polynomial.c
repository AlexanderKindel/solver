#include "declarations.h"

struct NestedPolynomial*nested_polynomial_copy(struct Stack*output_stack, struct NestedPolynomial*a)
{
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_polynomial_copy(output_stack, a->coefficients[i]);
    }
    return out;
}

void nested_polynomial_trim_leading_zeroes(struct NestedPolynomial*a)
{
    for (size_t i = a->coefficient_count; i-- > 0;)
    {
        if (rational_polynomial_equals(a->coefficients[i], polynomial_zero))
        {
            a->coefficient_count -= 1;
        }
        else
        {
            return;
        }
    }
}

bool nested_polynomial_equals(struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    if (a->coefficient_count != b->coefficient_count)
    {
        return false;
    }
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (!rational_polynomial_equals(a->coefficients[i], b->coefficients[i]))
        {
            return false;
        }
    }
    return true;
}

struct NestedPolynomial*nested_polynomial_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        SWAP(a, b, struct NestedPolynomial*);
    }
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_polynomial_add(output_stack, local_stack,
            a->coefficients[i], b->coefficients[i]);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_polynomial_copy(output_stack, a->coefficients[i]);
    }
    nested_polynomial_trim_leading_zeroes(out);
    return out;
}

struct NestedPolynomial*nested_polynomial_negate(struct Stack*restrict output_stack,
    struct NestedPolynomial*a)
{
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_polynomial_negate(output_stack, a->coefficients[i]);
    }
    return out;
}

struct NestedPolynomial*nested_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*minuend,
    struct NestedPolynomial*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*out = nested_polynomial_add(output_stack, local_stack, minuend,
        nested_polynomial_negate(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*nested_polynomial_rational_polynomial_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*b)
{
    if (rational_polynomial_equals(b, polynomial_zero))
    {
        return polynomial_zero;
    }
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            rational_polynomial_multiply(output_stack, local_stack, a->coefficients[i], b);
    }
    return out;
}

struct NestedPolynomial*nested_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_zero;
    }
    struct NestedPolynomial*out = POLYNOMIAL_ALLOCATE(output_stack,
        a->coefficient_count + b->coefficient_count - 1, struct RationalPolynomial*);
    for (size_t i = 0; i < out->coefficient_count; ++i)
    {
        out->coefficients[i] = polynomial_zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            out->coefficients[i + j] =
                rational_polynomial_add(output_stack, local_stack, out->coefficients[i + j],
                    rational_polynomial_multiply(local_stack, output_stack, a->coefficients[i],
                        b->coefficients[j]));
        }
    }
    nested_polynomial_trim_leading_zeroes(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*nested_polynomial_rational_polynomial_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*dividend, struct RationalPolynomial*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, dividend->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < dividend->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_polynomial_euclidean_divide(output_stack, local_stack,
            dividend->coefficients[i], divisor).quotient;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void nested_polynomial_copy_coefficients(struct Stack*output_stack, struct NestedPolynomial*a)
{
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        a->coefficients[i] = rational_polynomial_copy(output_stack, a->coefficients[i]);
    }
}

//Sets out.quotient and out.remainder to 0 when the coefficients of the quotient wouldn't all be
//rational polynomials.
struct NestedPolynomialDivision nested_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*dividend, struct NestedPolynomial*divisor)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        return (struct NestedPolynomialDivision) { polynomial_zero,
            nested_polynomial_copy(output_stack, dividend) };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomialDivision out = { POLYNOMIAL_ALLOCATE(output_stack,
        1 + dividend->coefficient_count - divisor->coefficient_count, struct RationalPolynomial*),
        POLYNOMIAL_ALLOCATE(output_stack, dividend->coefficient_count,
            struct RationalPolynomial*) };
    memcpy(out.remainder->coefficients, dividend->coefficients,
        dividend->coefficient_count * sizeof(struct RationalPolynomial*));
    while (out.remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out.remainder->coefficient_count;
        struct RationalPolynomialDivision division =
            rational_polynomial_euclidean_divide(local_stack, output_stack,
                out.remainder->coefficients[out.remainder->coefficient_count],
                divisor->coefficients[divisor->coefficient_count - 1]);
        if (rational_polynomial_equals(division.remainder, polynomial_zero))
        {
            out.quotient->coefficients[out.remainder->coefficient_count +
                1 - divisor->coefficient_count] = division.quotient;
            for (size_t i = 1; i < divisor->coefficient_count; ++i)
            {
                out.remainder->coefficients[out.remainder->coefficient_count - i] =
                    rational_polynomial_subtract(local_stack, output_stack,
                        out.remainder->coefficients[out.remainder->coefficient_count - i],
                        rational_polynomial_multiply(local_stack, output_stack, division.quotient,
                            divisor->coefficients[divisor->coefficient_count - i - 1]));
            }
        }
        else
        {
            out.quotient = 0;
            out.remainder = 0;
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
    }
    nested_polynomial_copy_coefficients(output_stack, out.quotient);
    nested_polynomial_trim_leading_zeroes(out.remainder);
    nested_polynomial_copy_coefficients(output_stack, out.remainder);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*nested_polynomial_get_content(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return rational_polynomial_one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = a->coefficients[0];
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        out = rational_polynomial_get_gcd(local_stack, output_stack, out, a->coefficients[i]);
    }
    out = rational_polynomial_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Correct up to multiplication by a nonzero rational.
struct RationalPolynomial*nested_polynomial_get_resultant(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b)
{
    if (a->coefficient_count == 0 || b->coefficient_count == 0)
    {
        return polynomial_zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*a_content =
        nested_polynomial_get_content(local_stack, output_stack, a);
    struct RationalPolynomial*b_content =
        nested_polynomial_get_content(local_stack, output_stack, b);
    struct RationalPolynomial*t = rational_polynomial_multiply(local_stack, output_stack,
        rational_polynomial_exponentiate(local_stack, output_stack, a_content,
            size_t_to_integer(local_stack, b->coefficient_count - 1)),
        rational_polynomial_exponentiate(local_stack, output_stack, b_content,
            size_t_to_integer(local_stack, a->coefficient_count - 1)));
    a = nested_polynomial_rational_polynomial_divide(local_stack, output_stack, a, a_content);
    b = nested_polynomial_rational_polynomial_divide(local_stack, output_stack, b, b_content);
    struct RationalPolynomial*g = rational_polynomial_one;
    struct RationalPolynomial*h = rational_polynomial_one;
    if (b->coefficient_count > a->coefficient_count)
    {
        SWAP(a, b, struct NestedPolynomial*);
    }
    while (b->coefficient_count > 1)
    {
        size_t degree_size_t = a->coefficient_count - b->coefficient_count;
        struct Integer*degree = size_t_to_integer(local_stack, degree_size_t);
        struct NestedPolynomial*remainder =
            nested_polynomial_euclidean_divide(local_stack, output_stack,
                nested_polynomial_rational_polynomial_multiply(local_stack, output_stack, a,
                    rational_polynomial_exponentiate(local_stack, output_stack,
                        b->coefficients[b->coefficient_count - 1],
                        integer_add(local_stack, degree, &one))), b).remainder;
        a = b;
        b = nested_polynomial_rational_polynomial_divide(local_stack, output_stack, remainder,
            rational_polynomial_multiply(local_stack, output_stack, g,
                rational_polynomial_exponentiate(local_stack, output_stack, h, degree)));
        g = a->coefficients[a->coefficient_count - 1];
        struct RationalPolynomial*term =
            rational_polynomial_exponentiate(local_stack, output_stack, g, degree);
        if (degree_size_t > 1)
        {
            h = rational_polynomial_euclidean_divide(local_stack, output_stack, term,
                rational_polynomial_exponentiate(local_stack, output_stack, h,
                    size_t_to_integer(local_stack, degree_size_t - 1))).quotient;
        }
        else
        {
            h = rational_polynomial_multiply(local_stack, output_stack, term,
                rational_polynomial_exponentiate(local_stack, output_stack, h,
                    size_t_to_integer(local_stack, 1 - degree_size_t)));
        }
    }
    struct RationalPolynomial*term = rational_polynomial_exponentiate(local_stack, output_stack,
        b->coefficients[b->coefficient_count - 1],
        size_t_to_integer(local_stack, a->coefficient_count - 1));
    struct RationalPolynomial*out;
    if (a->coefficient_count > 2)
    {
        out = rational_polynomial_euclidean_divide(local_stack, output_stack, term,
            rational_polynomial_exponentiate(local_stack, output_stack, h,
                size_t_to_integer(local_stack, a->coefficient_count - 2))).quotient;
    }
    else
    {
        out = rational_polynomial_multiply(local_stack, output_stack, term,
            rational_polynomial_exponentiate(local_stack, output_stack, h,
                size_t_to_integer(local_stack, 2 - a->coefficient_count)));
    }
    out = rational_polynomial_multiply(output_stack, local_stack, out, t);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*nested_polynomial_get_derivative(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count - 1, struct RationalPolynomial*);
    struct Integer*multiplier = &zero;
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        multiplier = integer_add(local_stack, multiplier, &one);
        out->coefficients[i - 1] = rational_polynomial_rational_multiply(output_stack, local_stack,
            a->coefficients[i], &(struct Rational) { multiplier, &one });
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}