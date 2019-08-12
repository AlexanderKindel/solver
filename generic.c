#include "declarations.h"

__declspec(noreturn) void crash(char*message)
{
    puts(message);
    abort();
}

void array_reverse(void**a, size_t element_count)
{
    for (size_t i = 0; i < element_count / 2; ++i)
    {
        a[i] = a[element_count - 1 - i];
    }
}

void*generic_exponentiate(struct RingOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, void*base, struct Integer*exponent, void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    void*out = operations->multiplicative_identity;
    while (exponent->sign > 0)
    {
        if (exponent->value[0] & 1)
        {
            out = operations->multiply(local_stack, output_stack, out, base, misc);
        }
        base = operations->multiply(local_stack, output_stack, base, base, misc);
        exponent = integer_half(local_stack, exponent);
    }
    out = operations->copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*generic_gcd(struct EuclideanDomainOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, void*a, void*b, void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!operations->ring_operations.equals(b, operations->ring_operations.additive_identity))
    {
        void*c = b;
        struct Division division;
        operations->euclidean_divide(local_stack, output_stack, &division, a, b, misc);
        b = division.remainder;
        a = c;
    }
    a = operations->ring_operations.copy(output_stack, a);
    local_stack->cursor = local_stack_savepoint;
    return a;
}

void generic_extended_gcd(struct EuclideanDomainOperations*operations, struct Stack*output_stack,
    struct Stack*local_stack, struct ExtendedGCDInfo*out, void*a, void*b, void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->a_coefficient = operations->ring_operations.additive_identity;
    out->b_coefficient = operations->ring_operations.multiplicative_identity;
    out->a_over_gcd = operations->ring_operations.additive_identity;
    out->b_over_gcd = operations->ring_operations.multiplicative_identity;
    while (!operations->ring_operations.equals(a, operations->ring_operations.additive_identity))
    {
        struct Division division;
        operations->euclidean_divide(local_stack, output_stack, &division, b, a, misc);
        void*m = operations->ring_operations.add(local_stack, output_stack, out->a_coefficient,
            operations->ring_operations.negative(local_stack, output_stack,
                operations->ring_operations.multiply(local_stack, output_stack, out->b_over_gcd,
                    division.quotient, misc), misc), misc);
        void*n = operations->ring_operations.add(local_stack, output_stack, out->b_coefficient,
            operations->ring_operations.negative(local_stack, output_stack,
                operations->ring_operations.multiply(local_stack, output_stack, out->a_over_gcd,
                    division.quotient, misc), misc), misc);
        b = a;
        a = division.remainder;
        out->a_coefficient = out->b_over_gcd;
        out->b_coefficient = out->a_over_gcd;
        out->b_over_gcd = m;
        out->a_over_gcd = n;
    }
    out->a_coefficient = operations->ring_operations.copy(output_stack, out->a_coefficient);
    out->a_over_gcd = operations->ring_operations.copy(output_stack, out->a_over_gcd);
    out->b_coefficient = operations->ring_operations.copy(output_stack, out->b_coefficient);
    out->b_over_gcd =
        operations->ring_operations.negative(output_stack, local_stack, out->b_over_gcd, misc);
    out->gcd = operations->ring_operations.copy(output_stack, b);
    local_stack->cursor = local_stack_savepoint;
}

size_t polynomial_size(size_t coefficient_count)
{
    return sizeof(struct Polynomial) + coefficient_count * sizeof(size_t);
}

void*polynomial_allocate(struct Stack*output_stack, size_t coefficient_count)
{
    struct Polynomial*out = stack_slot_allocate(output_stack,
        polynomial_size(coefficient_count), _Alignof(struct Polynomial));
    out->coefficient_count = coefficient_count;
    return out;
}

void polynomial_copy_coefficients(void*(coefficient_copy)(struct Stack*, void*),
    struct Stack*output_stack, struct Polynomial*a)
{
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        a->coefficients[i] = coefficient_copy(output_stack, a->coefficients[i]);
    }
}

void*polynomial_copy(void*(coefficient_copy)(struct Stack*, void*), struct Stack*output_stack,
    struct Polynomial*a)
{
    struct Polynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    memcpy(out->coefficients, a->coefficients, a->coefficient_count * sizeof(void*));
    polynomial_copy_coefficients(coefficient_copy, output_stack, out);
    return out;
}

void polynomial_trim_leading_zeroes(struct RingOperations*coefficient_operations,
    struct Polynomial*a)
{
    for (size_t i = a->coefficient_count; i-- > 0;)
    {
        if (coefficient_operations->equals(a->coefficients[i],
            coefficient_operations->additive_identity))
        {
            a->coefficient_count -= 1;
        }
        else
        {
            return;
        }
    }
}

bool polynomial_equals(struct RingOperations*coefficient_operations, struct Polynomial*a,
    struct Polynomial*b)
{
    if (a->coefficient_count != b->coefficient_count)
    {
        return false;
    }
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (!coefficient_operations->equals(a->coefficients[i], b->coefficients[i]))
        {
            return false;
        }
    }
    return true;
}

void*polynomial_add(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        POINTER_SWAP(a, b);
    }
    struct Polynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] = coefficient_operations->add(output_stack, local_stack,
            a->coefficients[i], b->coefficients[i], 0);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = coefficient_operations->copy(output_stack, a->coefficients[i]);
    }
    polynomial_trim_leading_zeroes(coefficient_operations, out);
    return out;
}

void*polynomial_negative(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Polynomial*a)
{
    struct Polynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            coefficient_operations->negative(output_stack, 0, a->coefficients[i], 0);
    }
    return out;
}

void*polynomial_subtract(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*minuend, struct Polynomial*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    void*out = polynomial_add(coefficient_operations, output_stack, local_stack, minuend,
        polynomial_negative(coefficient_operations, local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*polynomial_multiply(struct RingOperations*coefficient_operations, struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial*b, void*misc)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_allocate(output_stack, 0);
    }
    struct Polynomial*out =
        polynomial_allocate(output_stack, a->coefficient_count + b->coefficient_count - 1);
    for (size_t i = 0; i < out->coefficient_count; ++i)
    {
        out->coefficients[i] = coefficient_operations->additive_identity;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            out->coefficients[i + j] = coefficient_operations->add(output_stack, local_stack,
                out->coefficients[i + j], coefficient_operations->multiply(local_stack,
                    output_stack, a->coefficients[i], b->coefficients[j], misc), misc);
        }
    }
    polynomial_trim_leading_zeroes(coefficient_operations, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*polynomial_multiply_by_coefficient(struct RingOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a, void*b, void*misc)
{
    if (coefficient_operations->equals(b, coefficient_operations->additive_identity))
    {
        return &polynomial_zero;
    }
    struct Polynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = coefficient_operations->multiply(output_stack, local_stack,
            a->coefficients[i], b, misc);
    }
    return out;
}

//Sets out->quotient to 0 when the coefficients of the quotient wouldn't be in the ring.
void polynomial_euclidean_divide(struct EuclideanDomainOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct PolynomialDivision*out,
    struct Polynomial*dividend, struct Polynomial*divisor, void*misc)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        out->quotient = &polynomial_zero;
        out->remainder =
            polynomial_copy(coefficient_operations->ring_operations.copy, output_stack, dividend);
        return;
    }
    void*local_stack_savepoint = local_stack->cursor;
    out->quotient = polynomial_allocate(output_stack,
        1 + dividend->coefficient_count - divisor->coefficient_count);
    out->remainder = polynomial_allocate(output_stack, dividend->coefficient_count);
    memcpy(out->remainder->coefficients, dividend->coefficients,
        dividend->coefficient_count * sizeof(void*));
    while (out->remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out->remainder->coefficient_count;
        struct Division division;
        coefficient_operations->euclidean_divide(local_stack, output_stack, &division,
            out->remainder->coefficients[out->remainder->coefficient_count],
            divisor->coefficients[divisor->coefficient_count - 1], misc);
        if (coefficient_operations->ring_operations.equals(division.remainder,
            coefficient_operations->ring_operations.additive_identity))
        {
            out->quotient->coefficients[out->remainder->coefficient_count + 1 -
                divisor->coefficient_count] = division.quotient;
            for (size_t i = 1; i < divisor->coefficient_count; ++i)
            {
                out->remainder->coefficients[out->remainder->coefficient_count - i] =
                    coefficient_operations->ring_operations.add(local_stack, output_stack,
                        out->remainder->coefficients[out->remainder->coefficient_count - i],
                        coefficient_operations->ring_operations.negative(local_stack, output_stack,
                            coefficient_operations->ring_operations.multiply(local_stack,
                                output_stack, division.quotient,
                                divisor->coefficients[divisor->coefficient_count - i - 1], misc),
                            misc), misc);
            }
        }
        else
        {
            out->quotient = 0;
            local_stack->cursor = local_stack_savepoint;
            return;
        }
    }
    polynomial_copy_coefficients(coefficient_operations->ring_operations.copy, output_stack,
        out->quotient);
    polynomial_trim_leading_zeroes(&coefficient_operations->ring_operations, out->remainder);
    polynomial_copy_coefficients(coefficient_operations->ring_operations.copy, output_stack,
        out->remainder);
    local_stack->cursor = local_stack_savepoint;
}

void field_polynomial_euclidean_divide(struct FieldOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct PolynomialDivision*out,
    struct Polynomial*dividend, struct Polynomial*divisor, void*misc)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        out->quotient = polynomial_allocate(output_stack, 0);
        out->remainder =
            polynomial_copy(coefficient_operations->ring_operations.copy, output_stack, dividend);
        return;
    }
    void*local_stack_savepoint = local_stack->cursor;
    out->quotient = polynomial_allocate(output_stack,
        1 + dividend->coefficient_count - divisor->coefficient_count);
    out->remainder = polynomial_allocate(output_stack, dividend->coefficient_count);
    memcpy(out->remainder->coefficients, dividend->coefficients,
        dividend->coefficient_count * sizeof(void*));
    void*leading_coefficient_reciprocal = coefficient_operations->reciprocal(local_stack,
        output_stack, divisor->coefficients[divisor->coefficient_count - 1], misc);
    while (out->remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out->remainder->coefficient_count;
        void*quotient = coefficient_operations->ring_operations.multiply(local_stack, output_stack,
            out->remainder->coefficients[out->remainder->coefficient_count],
            leading_coefficient_reciprocal, misc);
        out->quotient->coefficients[out->remainder->coefficient_count + 1 -
            divisor->coefficient_count] = quotient;
        for (size_t i = 1; i < divisor->coefficient_count; ++i)
        {
            out->remainder->coefficients[out->remainder->coefficient_count - i] =
                coefficient_operations->ring_operations.add(local_stack, output_stack,
                    out->remainder->coefficients[out->remainder->coefficient_count - i],
                    coefficient_operations->ring_operations.negative(local_stack, output_stack,
                        coefficient_operations->ring_operations.multiply(local_stack, output_stack,
                            quotient, divisor->coefficients[divisor->coefficient_count - i - 1],
                            misc), misc), misc);
        }
    }
    polynomial_copy_coefficients(coefficient_operations->ring_operations.copy, output_stack,
        out->quotient);
    polynomial_trim_leading_zeroes(&coefficient_operations->ring_operations, out->remainder);
    polynomial_copy_coefficients(coefficient_operations->ring_operations.copy, output_stack,
        out->remainder);
    local_stack->cursor = local_stack_savepoint;
}

void*polynomial_divide_by_coefficient(void*(coefficient_quotient)(struct Stack*, struct Stack*,
    void*, void*), struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*dividend,
    void*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Polynomial*out = polynomial_allocate(output_stack, dividend->coefficient_count);
    for (size_t i = 0; i < dividend->coefficient_count; ++i)
    {
        out->coefficients[i] =
            coefficient_quotient(output_stack, local_stack, dividend->coefficients[i], divisor);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*polynomial_content(struct EuclideanDomainOperations*coefficient_operations,
    struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a, void*misc)
{
    if (!a->coefficient_count)
    {
        return coefficient_operations->ring_operations.multiplicative_identity;
    }
    void*local_stack_savepoint = local_stack->cursor;
    void*out = a->coefficients[0];
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        out = coefficient_operations->gcd(local_stack, output_stack, out, a->coefficients[i], misc);
    }
    out = coefficient_operations->ring_operations.copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t polynomial_squarefree_factor(struct EuclideanDomainOperations*polynomial_operations,
    void*(derivative)(struct Stack*, struct Stack*, void*), struct Stack*output_stack,
    struct Stack*local_stack, struct Polynomial*a, struct Polynomial**out, void*misc)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    struct Polynomial*b = a;
    struct Polynomial*c = derivative(local_stack, output_stack, a);
    a = polynomial_operations->gcd(local_stack, output_stack, b, c, misc);
    while (true)
    {
        struct Division division;
        polynomial_operations->euclidean_divide(local_stack, output_stack, &division, b, a, misc);
        b = division.quotient;
        if (polynomial_operations->ring_operations.equals(b,
            polynomial_operations->ring_operations.multiplicative_identity))
        {
            local_stack->cursor = local_stack_savepoint;
            return factor_count;
        }
        polynomial_operations->euclidean_divide(local_stack, output_stack, &division, c, a, misc);
        c = polynomial_operations->ring_operations.add(local_stack, output_stack, division.quotient,
            polynomial_operations->ring_operations.negative(local_stack, output_stack,
                derivative(local_stack, output_stack, b), misc), misc);
        a = polynomial_operations->gcd(output_stack, local_stack, b, c, misc);
        out[factor_count] = a;
        ++factor_count;
    }
}

void*polynomial_derivative(void*(coefficient_times_integer)(struct Stack*, struct Stack*, void*,
    struct Integer*), struct Stack*output_stack, struct Stack*local_stack, struct Polynomial*a)
{
    if (!a->coefficient_count)
    {
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Polynomial*out = polynomial_allocate(output_stack, a->coefficient_count - 1);
    struct Integer*multiplier = &zero;
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        multiplier = integer_add(local_stack, multiplier, &one);
        out->coefficients[i - 1] =
            coefficient_times_integer(output_stack, local_stack, a->coefficients[i], multiplier);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void*integer_generic_add(struct Stack*output_stack, struct Stack*unused_stack, void*a, void*b,
    void*unused)
{
    return integer_add(output_stack, a, b);
}

void*integer_generic_negative(struct Stack*output_stack, struct Stack*unused_stack, void*a,
    void*unused)
{
    return integer_negative(output_stack, a);
}

void*integer_generic_multiply(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused)
{
    return integer_multiply(output_stack, local_stack, a, b);
}

void integer_generic_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Division*out, void*dividend, void*divisor, void*unused)
{
    integer_euclidean_divide(output_stack, local_stack, (struct IntegerDivision*)out, dividend,
        divisor);
}

void*integer_generic_gcd(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused)
{
    return integer_gcd(output_stack, local_stack, a, b);
}

void*rational_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused)
{
    return rational_add(output_stack, local_stack, a, b);
}

void*rational_generic_negative(struct Stack*output_stack, struct Stack*unused_stack, void*a,
    void*unused)
{
    return rational_negative(output_stack, a);
}

void*rational_generic_multiply(struct Stack*output_stack, struct Stack*local_stack, void*a, void*b,
    void*unused)
{
    return rational_multiply(output_stack, local_stack, a, b);
}

void*rational_generic_reciprocal(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*unused)
{
    return rational_reciprocal(output_stack, a);
}

void*rational_generic_divide(struct Stack*output_stack, struct Stack*local_stack, void*dividend,
    void*divisor, void*unused)
{
    return rational_divide(output_stack, local_stack, dividend, divisor);
}

struct Float*float_generic_add(struct Stack*output_stack, struct Stack*local_stack, struct Float*a,
    struct Float*b, void*unused)
{
    return float_add(output_stack, local_stack, a, b);
}

struct Float*float_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    struct Float*a, void*unused)
{
    return float_negative(output_stack, a);
}

struct Float*float_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Float*a, struct Float*b, void*unused)
{
    return float_multiply(output_stack, local_stack, a, b);
}

void*gaussian_rational_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused)
{
    return gaussian_rational_add(output_stack, local_stack, a, b);
}

void*gaussian_rational_generic_multiply(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused)
{
    return gaussian_rational_multiply(output_stack, local_stack, a, b);
}

void*integer_polynomial_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused)
{
    return integer_polynomial_add(output_stack, local_stack, a, b);
}

void*integer_polynomial_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    void*a, void*unused)
{
    return integer_polynomial_negative(output_stack, a);
}

void*integer_polynomial_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    void*a, void*b, void*unused)
{
    return integer_polynomial_multiply(output_stack, local_stack, a, b);
}

void integer_polynomial_generic_euclidean_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct Division*out, void*dividend, void*divisor, void*unused)
{
    integer_polynomial_euclidean_divide(output_stack, local_stack, (struct PolynomialDivision*)out,
        dividend, divisor);
}

void*integer_polynomial_generic_gcd(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused)
{
    return integer_polynomial_gcd(output_stack, local_stack, a, b);
}

void*rational_polynomial_generic_add(struct Stack*output_stack, struct Stack*local_stack, void*a,
    void*b, void*unused)
{
    return rational_polynomial_add(output_stack, local_stack, a, b);
}

void*rational_polynomial_generic_negative(struct Stack*output_stack, struct Stack*unused_stack,
    void*a, void*unused)
{
    return rational_polynomial_negative(output_stack, a);
}

void*rational_polynomial_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    void*a, void*b, void*unused)
{
    return rational_polynomial_multiply(output_stack, local_stack, a, b);
}

void rational_polynomial_generic_euclidean_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct Division*out, void*dividend, void*divisor, void*unused)
{
    rational_polynomial_euclidean_divide(output_stack, local_stack, (struct PolynomialDivision*)out,
        dividend, divisor);
}

void*gaussian_rational_polynomial_generic_multiply(struct Stack*output_stack,
    struct Stack*local_stack, void*a, void*b, void*unused)
{
    return gaussian_rational_polynomial_multiply(output_stack, local_stack, a, b);
}

void*nested_polynomial_generic_add(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, struct NestedPolynomial*b, void*unused)
{
    return nested_polynomial_add(output_stack, local_stack, a, b);
}

void*nested_polynomial_generic_negative(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, void*unused)
{
    return nested_polynomial_negative(output_stack, local_stack, a);
}

void*number_generic_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Number*b, void*unused)
{
    return number_multiply(output_stack, local_stack, a, b);
}