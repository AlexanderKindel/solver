#include "declarations.h"

struct GaussianRational*gaussian_rational_copy(struct Stack*output_stack, struct GaussianRational*a)
{
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_copy(output_stack, a->real);
    out->imaginary = rational_copy(output_stack, a->imaginary);
    return out;
}

bool gaussian_rational_equals(struct GaussianRational*a, struct GaussianRational*b)
{
    return rational_equals(a->real, b->real) && rational_equals(a->imaginary, b->imaginary);
}

struct GaussianRational*gaussian_rational_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b)
{
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_add(output_stack, local_stack, a->real, b->real);
    out->imaginary = rational_add(output_stack, local_stack, a->imaginary, b->imaginary);
    return out;
}

struct GaussianRational*gaussian_rational_negate(struct Stack*output_stack,
    struct GaussianRational*a)
{
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_negate(output_stack, a->real);
    out->imaginary = rational_negate(output_stack, a->imaginary);
    return out;
}

struct GaussianRational*gaussian_rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*minuend,
    struct GaussianRational*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational*out = gaussian_rational_add(output_stack, local_stack, minuend,
        gaussian_rational_negate(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRational*gaussian_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_subtract(output_stack, local_stack,
        rational_multiply(local_stack, output_stack, a->real, b->real),
        rational_multiply(local_stack, output_stack, a->imaginary, b->imaginary));
    out->imaginary = rational_add(output_stack, local_stack,
        rational_multiply(local_stack, output_stack, a->real, b->imaginary),
        rational_multiply(local_stack, output_stack, a->imaginary, b->real));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRational*gaussian_rational_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct Rational*b)
{
    return gaussian_rational_multiply(output_stack, local_stack, a,
        &(struct GaussianRational){b, &rational_zero});
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_rational_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct GaussianRationalPolynomial*a, struct Rational*b)
{
    return gaussian_rational_polynomial_gaussian_rational_multiply(output_stack, local_stack, a,
        &(struct GaussianRational){b, &rational_zero});
}