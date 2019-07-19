#include "declarations.h"

struct Integer*modded_integer_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Integer*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo info;
    integer_extended_gcd(local_stack, output_stack, &info, a, characteristic);
    struct Integer*out;
    if (((struct Integer*)info.a_coefficient)->sign < 0)
    {
        out = integer_add(output_stack, info.a_coefficient, characteristic);
    }
    else
    {
        out = integer_copy_to_stack(output_stack, info.a_coefficient);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_reduced(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    struct IntegerPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_euclidean_remainder(local_stack, output_stack,
            a->coefficients[i], characteristic);
        if (out->coefficients[i]->sign < 0)
        {
            out->coefficients[i] = integer_add(local_stack, out->coefficients[i], characteristic);
        }
    }
    polynomial_trim_leading_zeroes(&integer_operations.ring_operations, (struct Polynomial*)out);
    polynomial_copy_coefficients(integer_copy_to_stack, output_stack, (struct Polynomial*)out);
    return out;
}

struct IntegerPolynomial*modded_polynomial_add(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_add(local_stack, output_stack, a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_negative(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_negative(local_stack, a), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_subtract(local_stack, output_stack, a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        integer_polynomial_multiply(local_stack, output_stack, a, b), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_multiply_by_coefficient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduced(output_stack, local_stack,
        polynomial_multiply_by_coefficient(&modded_integer_operations.ring_operations, local_stack,
            output_stack, (struct Polynomial*)a, b, 0), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void modded_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    field_polynomial_euclidean_divide(&modded_integer_operations, local_stack, output_stack, out,
        (struct Polynomial*)dividend, (struct Polynomial*)divisor, characteristic);
    out->quotient = (struct Polynomial*)modded_polynomial_reduced(output_stack, local_stack,
        (struct IntegerPolynomial*)out->quotient, characteristic);
    out->remainder = (struct Polynomial*)modded_polynomial_reduced(output_stack, local_stack,
        (struct IntegerPolynomial*)out->remainder, characteristic);
    local_stack->cursor = local_stack_savepoint;
}

struct IntegerPolynomial*modded_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic)
{
    struct PolynomialDivision division;
    modded_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor,
        characteristic);
    return (struct IntegerPolynomial*)division.quotient;
}

struct IntegerPolynomial*modded_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic)
{
    struct PolynomialDivision division;
    modded_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor,
        characteristic);
    return (struct IntegerPolynomial*)division.remainder;
}

struct IntegerPolynomial*modded_polynomial_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_multiply_by_coefficient(output_stack,
        local_stack, a, modded_integer_reciprocal(local_stack, output_stack,
            a->coefficients[a->coefficient_count - 1], characteristic), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_exponentiate(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*base, struct Integer*exponent,
    struct Integer*characteristic)
{
    return generic_exponentiate(&modded_polynomial_operations.ring_operations, output_stack,
        local_stack, base, exponent, characteristic);
}

struct IntegerPolynomial*modded_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b, struct Integer*characteristic)
{
    return generic_gcd(&modded_polynomial_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b, characteristic);
}

void modded_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct ExtendedGCDInfo*out, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    generic_extended_gcd(&modded_polynomial_operations, output_stack, local_stack, out,
        (struct Polynomial*)a, (struct Polynomial*)b, characteristic);
}

size_t cantor_zassenhaus_split(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, size_t degree,
    struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    if ((a->coefficient_count - 1) / degree == 1)
    {
        out[0] = modded_polynomial_multiply_by_coefficient(output_stack, local_stack, a,
            modded_integer_reciprocal(local_stack, output_stack,
            a->coefficients[a->coefficient_count - 1], characteristic), characteristic);
        local_stack->cursor = local_stack_savepoint;
        return 1;
    }
    struct IntegerPolynomial*b;
    if (integer_equals(characteristic, &INT(2, +)))
    {
        struct IntegerPolynomial*x_squared = polynomial_allocate(local_stack, 3);
        x_squared->coefficients[0] = &zero;
        x_squared->coefficients[1] = &zero;
        x_squared->coefficients[2] = &one;
        struct IntegerPolynomial*t = &integer_polynomial_one;
        struct IntegerPolynomial*c = &integer_polynomial_one;
        while (true)
        {
            for (size_t i = 1; i < degree; ++i)
            {
                c = modded_polynomial_add(local_stack, output_stack, t,
                    modded_polynomial_euclidean_remainder(local_stack, output_stack,
                        modded_polynomial_multiply(local_stack, output_stack, c, c, characteristic),
                        a, characteristic), characteristic);
            }
            b = modded_polynomial_gcd(local_stack, output_stack, a, c, characteristic);
            if (b->coefficient_count > 1 && b->coefficient_count != a->coefficient_count)
            {
                break;
            }
            t = modded_polynomial_multiply(local_stack, output_stack, t, x_squared, characteristic);
            c = t;
        }
    }
    else
    {
        struct Integer*p_minus_one = integer_add(local_stack, characteristic, &INT(1, -));
        struct IntegerPolynomial*t = polynomial_allocate(local_stack, 2 * degree);
        for (size_t i = 0; i < t->coefficient_count - 1; ++i)
        {
            t->coefficients[i] = p_minus_one;
        }
        t->coefficients[t->coefficient_count - 1] = &one;
        struct Integer*power = integer_add(local_stack,
            integer_exponentiate(local_stack, output_stack, characteristic,
                integer_from_size_t(local_stack, degree)),
            &INT(1, -));
        integer_halve(power);
        while (true)
        {
            b = modded_polynomial_gcd(local_stack, output_stack, a,
                modded_polynomial_subtract(local_stack, output_stack,
                    modded_polynomial_exponentiate(local_stack, output_stack, t, power,
                        characteristic), &integer_polynomial_one, characteristic), characteristic);
            if (b->coefficient_count > 1 && b->coefficient_count != a->coefficient_count)
            {
                break;
            }
        }
        size_t i = t->coefficient_count - 2;
        while (true)
        {
            if (t->coefficients[i]->value_count)
            {
                t->coefficients[i] = integer_add(local_stack, t->coefficients[i], &INT(1, -));
                break;
            }
            if (i == 0)
            {
                t->coefficient_count -= 1;
                for (size_t i = 0; i < t->coefficient_count - 1; ++i)
                {
                    t->coefficients[i] = p_minus_one;
                }
                t->coefficients[t->coefficient_count - 1] = &one;
                break;
            }
            t->coefficients[i] = p_minus_one;
            --i;
        }
    }
    size_t factor_count =
        cantor_zassenhaus_split(output_stack, local_stack, b, characteristic, degree, out);
    factor_count += cantor_zassenhaus_split(output_stack, local_stack,
        modded_polynomial_euclidean_quotient(local_stack, output_stack, a, b, characteristic),
        characteristic, degree, out + factor_count);
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

size_t squarefree_modded_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    struct IntegerPolynomial*x = polynomial_allocate(local_stack, 2);
    x->coefficients[0] = &zero;
    x->coefficients[1] = &one;
    struct IntegerPolynomial*v = a;
    struct IntegerPolynomial*w = x;
    size_t d = 0;
    while (2 * d + 2 < a->coefficient_count)
    {
        ++d;
        w = modded_polynomial_euclidean_remainder(local_stack, output_stack,
            modded_polynomial_exponentiate(local_stack, output_stack, w, characteristic,
                characteristic),
            a, characteristic);
        struct IntegerPolynomial*degree_d_factor_product =
            modded_polynomial_gcd(local_stack, output_stack,
                modded_polynomial_subtract(local_stack, output_stack, w, x, characteristic),
                v, characteristic);
        if (degree_d_factor_product->coefficient_count > 1)
        {
            factor_count += cantor_zassenhaus_split(output_stack, local_stack,
                degree_d_factor_product, characteristic, d, out + factor_count);
            v = modded_polynomial_euclidean_quotient(local_stack, output_stack, v,
                degree_d_factor_product, characteristic);
            w = modded_polynomial_euclidean_remainder(local_stack, output_stack, w, v,
                characteristic);
        }
    }
    if (v->coefficient_count > 1)
    {
        factor_count += cantor_zassenhaus_split(output_stack, local_stack, v, characteristic,
            v->coefficient_count - 1, out + factor_count);
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}