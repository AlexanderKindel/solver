#include "declarations.h"

struct GaussianRationalPolynomial*gaussian_rational_polynomial_copy(struct Stack*output_stack,
    struct GaussianRationalPolynomial*a)
{
    struct GaussianRationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct GaussianRational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = gaussian_rational_copy(output_stack, a->coefficients + i);
    }
    return out;
}

bool gaussian_rational_polynomial_equals(struct GaussianRationalPolynomial*a,
    struct GaussianRationalPolynomial*b)
{
    if (a->coefficient_count != b->coefficient_count)
    {
        return false;
    }
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (!gaussian_rational_equals(a->coefficients + i, b->coefficients + i))
        {
            return false;
        }
    }
    return true;
}

void gaussian_rational_polynomial_trim_leading_zeroes(struct GaussianRationalPolynomial*a)
{
    for (size_t i = a->coefficient_count; i-- > 0;)
    {
        if (gaussian_rational_equals(a->coefficients + i, &gaussian_rational_zero))
        {
            a->coefficient_count -= 1;
        }
        else
        {
            return;
        }
    }
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_add(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct GaussianRationalPolynomial*b)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        SWAP(a, b, struct GaussianRationalPolynomial*);
    }
    struct GaussianRationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct GaussianRational);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] = gaussian_rational_add(output_stack, local_stack, a->coefficients + i,
            b->coefficients + i);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = gaussian_rational_copy(output_stack, a->coefficients + i);
    }
    gaussian_rational_polynomial_trim_leading_zeroes(out);
    return out;
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_rational_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct Rational*b)
{
    return gaussian_rational_polynomial_gaussian_rational_multiply(output_stack, local_stack, a,
        &(struct GaussianRational) { *b, rational_zero });
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_gaussian_rational_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct GaussianRational*b)
{
    if (gaussian_rational_equals(b, &gaussian_rational_zero))
    {
        return polynomial_zero;
    }
    struct GaussianRationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct GaussianRational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            gaussian_rational_multiply(output_stack, local_stack, a->coefficients + i, b);
    }
    return out;
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct GaussianRationalPolynomial*b)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_zero;
    }
    struct GaussianRationalPolynomial*out = POLYNOMIAL_ALLOCATE(output_stack,
        a->coefficient_count + b->coefficient_count - 1, struct GaussianRational);
    for (size_t i = 0; i < out->coefficient_count; ++i)
    {
        out->coefficients[i] = gaussian_rational_zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            struct GaussianRational product = gaussian_rational_multiply(local_stack, output_stack,
                a->coefficients + i, b->coefficients + j);
            out->coefficients[i + j] = gaussian_rational_add(output_stack, local_stack,
                out->coefficients + i + j, &product);
        }
    }
    gaussian_rational_polynomial_trim_leading_zeroes(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}