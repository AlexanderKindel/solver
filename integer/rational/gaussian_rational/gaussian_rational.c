#include "declarations.h"

struct GaussianRational gaussian_rational_copy(struct Stack*output_stack, struct GaussianRational*a)
{
    return (struct GaussianRational) { rational_copy(output_stack, &a->real),
        rational_copy(output_stack, &a->imaginary) };
}

bool gaussian_rational_equals(struct GaussianRational*a, struct GaussianRational*b)
{
    return rational_equals(&a->real, &b->real) && rational_equals(&a->imaginary, &b->imaginary);
}

struct GaussianRational gaussian_rational_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b)
{
    return (struct GaussianRational) { rational_add(output_stack, local_stack, &a->real, &b->real),
        rational_add(output_stack, local_stack, &a->imaginary, &b->imaginary) };
}

struct GaussianRational gaussian_rational_negate(struct Stack*output_stack,
    struct GaussianRational*a)
{
    return (struct GaussianRational) { rational_negate(output_stack, &a->real),
        rational_negate(output_stack, &a->imaginary) };
}

struct GaussianRational gaussian_rational_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*minuend,
    struct GaussianRational*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational negative_subtrahend = gaussian_rational_negate(local_stack, subtrahend);
    struct GaussianRational out =
        gaussian_rational_add(output_stack, local_stack, minuend, &negative_subtrahend);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRational gaussian_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct GaussianRational*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct GaussianRational out;
    struct Rational product_a = rational_multiply(local_stack, output_stack, &a->real, &b->real);
    struct Rational product_b =
        rational_multiply(local_stack, output_stack, &a->imaginary, &b->imaginary);
    out.real = rational_subtract(output_stack, local_stack, &product_a, &product_b);
    product_a = rational_multiply(local_stack, output_stack, &a->real, &b->imaginary);
    product_b = rational_multiply(local_stack, output_stack, &a->imaginary, &b->real);
    out.imaginary = rational_add(output_stack, local_stack, &product_a, &product_b);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct GaussianRational gaussian_rational_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct GaussianRational*a, struct Rational*b)
{
    return (struct GaussianRational) { rational_multiply(output_stack, local_stack, &a->real, b),
        rational_multiply(output_stack, local_stack, &a->imaginary, b) };
}

#include "integer/rational/gaussian_rational/gaussian_rational_polynomial/gaussian_rational_polynomial.c"