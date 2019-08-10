#include "declarations.h"

struct GaussianRational*gaussian_rational_copy(struct Stack*output_stack, struct GaussianRational*a)
{
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_copy(output_stack, a->real);
    out->imaginary = rational_copy(output_stack, a->imaginary);
    return out;
}

struct GaussianRational*gaussian_rational_add(struct Stack*output_stack, struct Stack*local_stack,
    struct GaussianRational*a, struct GaussianRational*b)
{
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_add(output_stack, local_stack, a->real, b->real);
    out->imaginary = rational_add(output_stack, local_stack, a->imaginary, b->imaginary);
    return out;
}

struct GaussianRational*gaussian_rational_negative(struct Stack*output_stack,
    struct GaussianRational*a)
{
    struct GaussianRational*out = ALLOCATE(output_stack, struct GaussianRational);
    out->real = rational_negative(output_stack, a->real);
    out->imaginary = rational_negative(output_stack, a->imaginary);
    return out;
}

struct GaussianRational*gaussian_rational_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRational*a, struct GaussianRational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational*out = gaussian_rational_add(output_stack, local_stack, a,
        gaussian_rational_negative(local_stack, b));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRational*gaussian_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRational*a, struct GaussianRational*b)
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

struct GaussianRationalPolynomial*gaussian_rational_polynomial_copy(struct Stack*output_stack,
    struct GaussianRationalPolynomial*a)
{
    return polynomial_copy(gaussian_rational_copy, output_stack, (struct Polynomial*)a);
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_add(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRationalPolynomial*a,
    struct GaussianRationalPolynomial*b)
{
    return polynomial_add(&gaussian_rational_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct GaussianRationalPolynomial*a,
    struct GaussianRationalPolynomial*b)
{
    return polynomial_multiply(&gaussian_rational_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b, 0);
}

struct GaussianRationalPolynomial*gaussian_rational_polynomial_rational_multiply(
    struct Stack*output_stack, struct Stack*local_stack, struct GaussianRationalPolynomial*a,
    struct Rational*b)
{
    return polynomial_multiply_by_coefficient(&gaussian_rational_operations, output_stack,
        local_stack, (struct Polynomial*)a, &(struct GaussianRational){b, &rational_zero}, 0);
}