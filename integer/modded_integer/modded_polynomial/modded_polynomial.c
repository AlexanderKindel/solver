#include "declarations.h"

struct IntegerPolynomial*modded_polynomial_reduce(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            integer_get_residue(local_stack, output_stack, a->coefficients[i], characteristic);
    }
    integer_polynomial_trim_leading_zeroes(out);
    integer_polynomial_copy_coefficients(output_stack, out);
    return out;
}

struct IntegerPolynomial*modded_polynomial_add(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        SWAP(a, b, struct IntegerPolynomial*);
    }
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] = modded_integer_add(output_stack, local_stack, a->coefficients[i],
            b->coefficients[i], characteristic);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_copy(output_stack, a->coefficients[i]);
    }
    integer_polynomial_trim_leading_zeroes(out);
    return out;
}

struct IntegerPolynomial*modded_polynomial_negate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = modded_integer_negate(output_stack, local_stack, a->coefficients[i],
            characteristic);
    }
    return out;
}

struct IntegerPolynomial*modded_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*minuend,
    struct IntegerPolynomial*subtrahend, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_reduce(output_stack, local_stack,
        integer_polynomial_subtract(local_stack, output_stack, minuend, subtrahend), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_modded_integer_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*a, struct Integer*b, struct Integer*characteristic)
{
    if (integer_equals(b, &zero))
    {
        return polynomial_zero;
    }
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = modded_integer_multiply(output_stack, local_stack,
            a->coefficients[i], b, characteristic);
    }
    return out;
}

struct IntegerPolynomial*modded_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_zero;
    }
    struct IntegerPolynomial*out = POLYNOMIAL_ALLOCATE(output_stack,
        a->coefficient_count + b->coefficient_count - 1, struct Integer*);
    for (size_t i = 0; i < out->coefficient_count; ++i)
    {
        out->coefficients[i] = &zero;
    }
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        for (size_t j = 0; j < b->coefficient_count; ++j)
        {
            out->coefficients[i + j] =
                modded_integer_add(output_stack, local_stack, out->coefficients[i + j],
                    modded_integer_multiply(local_stack, output_stack, a->coefficients[i],
                        b->coefficients[j], characteristic), characteristic);
        }
    }
    integer_polynomial_trim_leading_zeroes(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomialDivision modded_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor,
    struct Integer*characteristic)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        return (struct IntegerPolynomialDivision) { polynomial_zero,
            integer_polynomial_copy(output_stack, dividend) };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomialDivision out = { POLYNOMIAL_ALLOCATE(output_stack,
        1 + dividend->coefficient_count - divisor->coefficient_count, struct Integer*),
        POLYNOMIAL_ALLOCATE(output_stack, dividend->coefficient_count, struct Integer*) };
    memcpy(out.remainder->coefficients, dividend->coefficients,
        dividend->coefficient_count * sizeof(void*));
    struct Integer*leading_coefficient_reciprocal = modded_integer_get_reciprocal(local_stack,
        output_stack, divisor->coefficients[divisor->coefficient_count - 1], characteristic);
    while (out.remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out.remainder->coefficient_count;
        struct Integer*quotient = modded_integer_multiply(local_stack, output_stack,
            out.remainder->coefficients[out.remainder->coefficient_count],
            leading_coefficient_reciprocal, characteristic);
        out.quotient->coefficients[out.remainder->coefficient_count +
            1 - divisor->coefficient_count] = quotient;
        for (size_t i = 1; i < divisor->coefficient_count; ++i)
        {
            out.remainder->coefficients[out.remainder->coefficient_count - i] =
                modded_integer_subtract(local_stack, output_stack,
                    out.remainder->coefficients[out.remainder->coefficient_count - i],
                    modded_integer_multiply(local_stack, output_stack, quotient,
                        divisor->coefficients[divisor->coefficient_count - i - 1], characteristic),
                    characteristic);
        }
    }
    integer_polynomial_copy_coefficients(output_stack, out.quotient);
    integer_polynomial_trim_leading_zeroes(out.remainder);
    integer_polynomial_copy_coefficients(output_stack, out.remainder);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*base, struct Integer*exponent,
    struct Integer*characteristic)
{
    if (!exponent->value_count)
    {
        return integer_polynomial_one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_one;
    while (true)
    {
        if (exponent->value[0] & 1)
        {
            out = modded_polynomial_multiply(local_stack, output_stack, out, base, characteristic);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = integer_polynomial_copy(output_stack, out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = modded_polynomial_multiply(local_stack, output_stack, base, base, characteristic);
    }
}

//Leaves a significant amount of excess allocations on output_stack.
struct ExtendedGCDInfo modded_polynomial_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo out =
    { b, polynomial_zero, integer_polynomial_one, polynomial_zero, integer_polynomial_one };
    while (!integer_polynomial_equals(a, polynomial_zero))
    {
        struct IntegerPolynomialDivision division = modded_polynomial_euclidean_divide(output_stack,
            local_stack, out.gcd, a, characteristic);
        struct IntegerPolynomial*m =
            modded_polynomial_subtract(output_stack, local_stack, out.a_coefficient,
                modded_polynomial_multiply(local_stack, output_stack, out.b_over_gcd,
                    division.quotient, characteristic), characteristic);
        struct IntegerPolynomial*n =
            modded_polynomial_subtract(output_stack, local_stack, out.b_coefficient,
                modded_polynomial_multiply(local_stack, output_stack, out.a_over_gcd,
                    division.quotient, characteristic), characteristic);
        out.gcd = a;
        a = division.remainder;
        out.a_coefficient = out.b_over_gcd;
        out.b_coefficient = out.a_over_gcd;
        out.b_over_gcd = m;
        out.a_over_gcd = n;
    }
    out.b_over_gcd =
        modded_polynomial_negate(output_stack, local_stack, out.b_over_gcd, characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*modded_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!integer_polynomial_equals(b, polynomial_zero))
    {
        struct IntegerPolynomial*c = b;
        b = modded_polynomial_euclidean_divide(local_stack, output_stack, a, b,
            characteristic).remainder;
        a = c;
    }
    a = integer_polynomial_copy(output_stack, a);
    local_stack->cursor = local_stack_savepoint;
    return a;
}

struct IntegerPolynomial*modded_polynomial_get_monic(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = modded_polynomial_modded_integer_multiply(output_stack,
        local_stack, a, modded_integer_get_reciprocal(local_stack, output_stack,
            a->coefficients[a->coefficient_count - 1], characteristic), characteristic);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t cantor_zassenhaus_split(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*characteristic,
    size_t degree, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    if ((a->coefficient_count - 1) / degree == 1)
    {
        out[0] = modded_polynomial_modded_integer_multiply(output_stack, local_stack, a,
            modded_integer_get_reciprocal(local_stack, output_stack,
            a->coefficients[a->coefficient_count - 1], characteristic), characteristic);
        local_stack->cursor = local_stack_savepoint;
        return 1;
    }
    struct IntegerPolynomial*b;
    if (integer_equals(characteristic, INT(2, 1)))
    {
        struct IntegerPolynomial*x_squared = POLYNOMIAL_ALLOCATE(local_stack, 3, struct Integer*);
        x_squared->coefficients[0] = &zero;
        x_squared->coefficients[1] = &zero;
        x_squared->coefficients[2] = &one;
        struct IntegerPolynomial*t = integer_polynomial_one;
        struct IntegerPolynomial*c = integer_polynomial_one;
        while (true)
        {
            for (size_t i = 1; i < degree; ++i)
            {
                c = modded_polynomial_add(local_stack, output_stack, t,
                    modded_polynomial_euclidean_divide(local_stack, output_stack,
                        modded_polynomial_multiply(local_stack, output_stack, c, c, characteristic),
                        a, characteristic).remainder,
                    characteristic);
            }
            b = modded_polynomial_get_gcd(local_stack, output_stack, a, c, characteristic);
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
        struct Integer*p_minus_one = integer_add(local_stack, characteristic, INT(1, -1));
        struct IntegerPolynomial*t = POLYNOMIAL_ALLOCATE(local_stack, 2 * degree, struct Integer*);
        for (size_t i = 0; i < t->coefficient_count - 1; ++i)
        {
            t->coefficients[i] = p_minus_one;
        }
        t->coefficients[t->coefficient_count - 1] = &one;
        struct Integer*power = integer_add(local_stack,
            integer_exponentiate(local_stack, output_stack, characteristic,
                size_t_to_integer(local_stack, degree)),
            INT(1, -1));
        integer_halve_in_place(power);
        while (true)
        {
            b = modded_polynomial_get_gcd(local_stack, output_stack, a,
                modded_polynomial_subtract(local_stack, output_stack,
                    modded_polynomial_exponentiate(local_stack, output_stack, t, power,
                        characteristic), integer_polynomial_one, characteristic), characteristic);
            if (b->coefficient_count > 1 && b->coefficient_count != a->coefficient_count)
            {
                break;
            }
            size_t i = t->coefficient_count - 2;
            while (true)
            {
                if (t->coefficients[i]->value_count)
                {
                    t->coefficients[i] = integer_add(local_stack, t->coefficients[i], INT(1, -1));
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
    }
    size_t factor_count =
        cantor_zassenhaus_split(output_stack, local_stack, b, characteristic, degree, out);
    factor_count += cantor_zassenhaus_split(output_stack, local_stack,
        modded_polynomial_euclidean_divide(local_stack, output_stack, a, b,
            characteristic).quotient,
        characteristic, degree, out + factor_count);
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

size_t squarefree_modded_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial**out, struct IntegerPolynomial*a,
    struct Integer*characteristic)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    POLY(x, struct IntegerPolynomial, struct Integer*, 2, &zero, &one);
    struct IntegerPolynomial*v = a;
    struct IntegerPolynomial*w = &x.p;
    size_t d = 0;
    while (2 * d + 2 < a->coefficient_count)
    {
        ++d;
        w = modded_polynomial_euclidean_divide(local_stack, output_stack,
            modded_polynomial_exponentiate(local_stack, output_stack, w, characteristic,
                characteristic),
            a, characteristic).remainder;
        struct IntegerPolynomial*degree_d_factor_product =
            modded_polynomial_get_gcd(local_stack, output_stack,
                modded_polynomial_subtract(local_stack, output_stack, w, &x.p, characteristic),
                v, characteristic);
        if (degree_d_factor_product->coefficient_count > 1)
        {
            factor_count += cantor_zassenhaus_split(output_stack, local_stack,
                degree_d_factor_product, characteristic, d, out + factor_count);
            v = modded_polynomial_euclidean_divide(local_stack, output_stack, v,
                degree_d_factor_product, characteristic).quotient;
            w = modded_polynomial_euclidean_divide(local_stack, output_stack, w, v,
                characteristic).remainder;
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
