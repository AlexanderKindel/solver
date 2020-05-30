#include "declarations.h"

struct NestedPolynomial*number_field_polynomial_number_field_element_multiply(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    if (rational_polynomial_equals(b, polynomial_zero))
    {
        return polynomial_zero;
    }
    struct NestedPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct RationalPolynomial*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = number_field_element_multiply(output_stack, local_stack,
            a->coefficients[i], b, generator_minimal_polynomial);
    }
    return out;
}

struct NestedPolynomial*number_field_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    if (!a->coefficient_count && !b->coefficient_count)
    {
        return polynomial_zero;
    }
    struct NestedPolynomial*out = POLYNOMIAL_ALLOCATE(output_stack,
        a->coefficient_count + b->coefficient_count - 1, struct RationalPolynoimial*);
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
                    number_field_element_multiply(local_stack, output_stack, a->coefficients[i],
                        b->coefficients[j], generator_minimal_polynomial));
        }
    }
    nested_polynomial_trim_leading_zeroes(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomialDivision number_field_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct NestedPolynomial*dividend, struct NestedPolynomial*divisor,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    if (divisor->coefficient_count > dividend->coefficient_count)
    {
        return (struct NestedPolynomialDivision) {
            polynomial_zero,
                nested_polynomial_copy(output_stack, dividend)
        };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomialDivision out = { POLYNOMIAL_ALLOCATE(output_stack,
        1 + dividend->coefficient_count - divisor->coefficient_count, struct RationalPolynomial*),
        POLYNOMIAL_ALLOCATE(output_stack, dividend->coefficient_count,
            struct RationalPolynomial*) };
    memcpy(out.remainder->coefficients, dividend->coefficients,
        dividend->coefficient_count * sizeof(void*));
    struct RationalPolynomial*leading_coefficient_reciprocal =
        number_field_element_get_reciprocal(local_stack, output_stack,
            divisor->coefficients[divisor->coefficient_count - 1], generator_minimal_polynomial);
    while (out.remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out.remainder->coefficient_count;
        struct RationalPolynomial*quotient = number_field_element_multiply(local_stack,
            output_stack, out.remainder->coefficients[out.remainder->coefficient_count],
            leading_coefficient_reciprocal, generator_minimal_polynomial);
        out.quotient->coefficients[out.remainder->coefficient_count +
            1 - divisor->coefficient_count] = quotient;
        for (size_t i = 1; i < divisor->coefficient_count; ++i)
        {
            out.remainder->coefficients[out.remainder->coefficient_count - i] =
                rational_polynomial_subtract(local_stack, output_stack,
                    out.remainder->coefficients[out.remainder->coefficient_count - i],
                    number_field_element_multiply(local_stack, output_stack, quotient,
                        divisor->coefficients[divisor->coefficient_count - i - 1],
                        generator_minimal_polynomial));
        }
    }
    nested_polynomial_copy_coefficients(output_stack, out.quotient);
    nested_polynomial_trim_leading_zeroes(out.remainder);
    nested_polynomial_copy_coefficients(output_stack, out.remainder);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct NestedPolynomial*number_field_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    while (!nested_polynomial_equals(b, polynomial_zero))
    {
        struct NestedPolynomial*c = b;
        b = number_field_polynomial_euclidean_divide(local_stack, output_stack, a, b,
            generator_minimal_polynomial).remainder;
        a = c;
    }
    a = nested_polynomial_copy(output_stack, a);
    local_stack->cursor = local_stack_savepoint;
    return a;
}

size_t number_field_polynomial_squarefree_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a, struct NestedPolynomial**out,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    struct NestedPolynomial*b = a;
    struct NestedPolynomial*c = nested_polynomial_get_derivative(local_stack, output_stack, a);
    a = number_field_polynomial_get_gcd(local_stack, output_stack, b, c,
        generator_minimal_polynomial);
    struct NestedPolynomialDivision division = number_field_polynomial_euclidean_divide(local_stack,
        output_stack, b, a, generator_minimal_polynomial);
    b = division.quotient;
    do
    {
        division = number_field_polynomial_euclidean_divide(local_stack, output_stack, c, a,
            generator_minimal_polynomial);
        c = nested_polynomial_subtract(local_stack, output_stack, division.quotient,
            nested_polynomial_get_derivative(local_stack, output_stack, b));
        a = number_field_polynomial_get_gcd(output_stack, local_stack, b, c,
            generator_minimal_polynomial);
        if (a->coefficient_count > 1)
        {
            out[factor_count] = a;
            ++factor_count;
        }
        division = number_field_polynomial_euclidean_divide(local_stack, output_stack, b, a,
            generator_minimal_polynomial);
        b = division.quotient;
    } while (b->coefficient_count > 1);
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

size_t number_field_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial**out, struct NestedPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial**squarefree_factors =
        ARRAY_ALLOCATE(local_stack, a->coefficient_count - 1, struct NestedPolynomial*);
    size_t squarefree_factor_count = number_field_polynomial_squarefree_factor(local_stack,
        output_stack, a, squarefree_factors, generator_minimal_polynomial);
    struct NestedPolynomial*nested_minimal_polynomial =
        rational_polynomial_to_nested_polynomial(local_stack, generator_minimal_polynomial);
    size_t factor_count = 0;
    for (size_t i = 0; i < squarefree_factor_count; ++i)
    {
        struct Integer*k = &zero;
        while (true)
        {
            struct NestedPolynomial*power = nested_polynomial_one;
            struct NestedPolynomial*d =
                POLYNOMIAL_ALLOCATE(local_stack, 2, struct RationalPolynomial*);
            d->coefficients[0] = POLYNOMIAL_ALLOCATE(local_stack, 2, struct Rational);
            d->coefficients[0]->coefficients[0] = rational_zero;
            d->coefficients[0]->coefficients[1] = rational_one;
            d->coefficients[1] = POLYNOMIAL_ALLOCATE(local_stack, 1, struct Rational);
            d->coefficients[1]->coefficients[0] = (struct Rational) { k, &one };
            struct NestedPolynomial*e = polynomial_zero;
            for (size_t j = 0; j < squarefree_factors[i]->coefficient_count; ++j)
            {
                e = nested_polynomial_add(local_stack, output_stack, e,
                    nested_polynomial_multiply(local_stack, output_stack, power,
                        rational_polynomial_to_nested_polynomial(local_stack,
                            squarefree_factors[i]->coefficients[j])));
                power = nested_polynomial_multiply(local_stack, output_stack, power, d);
            }
            struct RationalPolynomial*resultant = nested_polynomial_get_resultant(local_stack,
                output_stack, nested_minimal_polynomial, e);
            if (rational_polynomial_get_gcd(local_stack, output_stack,
                rational_polynomial_get_derivative(local_stack, output_stack, resultant),
                resultant)->coefficient_count < 2)
            {
                struct IntegerPolynomial**resultant_factors = ARRAY_ALLOCATE(local_stack,
                    resultant->coefficient_count - 1, struct IntegerPolynomial*);
                size_t resultant_factor_count =
                    primitive_integer_polynomial_factor(local_stack, output_stack,
                        resultant_factors,
                        rational_polynomial_get_primitive_part(local_stack, output_stack,
                            resultant));
                k->sign *= -1;
                d = POLYNOMIAL_ALLOCATE(local_stack, 2, struct RationalPolynomial*);
                d->coefficients[0] = POLYNOMIAL_ALLOCATE(local_stack, 2, struct Rational);
                d->coefficients[0]->coefficients[0] = rational_zero;
                d->coefficients[0]->coefficients[1] = (struct Rational) { k, &one };
                d->coefficients[1] = rational_polynomial_one;
                for (size_t j = 0; j < resultant_factor_count; ++j)
                {
                    power = nested_polynomial_one;
                    e = polynomial_zero;
                    for (size_t k = 0; k < resultant_factors[j]->coefficient_count; ++k)
                    {
                        POLY(rational_polynomial, struct RationalPolynomial, struct Rational, 1,
                            { resultant_factors[j]->coefficients[k], &one });
                        e = nested_polynomial_add(local_stack, output_stack, e,
                            nested_polynomial_rational_polynomial_multiply(local_stack,
                                output_stack, power, &rational_polynomial.p));
                        power = number_field_polynomial_multiply(local_stack, output_stack, power,
                            d, generator_minimal_polynomial);
                    }
                    struct NestedPolynomial*factor = number_field_polynomial_get_gcd(local_stack,
                        output_stack, squarefree_factors[i], e, generator_minimal_polynomial);
                    out[factor_count] =
                        number_field_polynomial_number_field_element_multiply(output_stack, local_stack,
                            factor,
                            number_field_element_get_reciprocal(local_stack, output_stack,
                                factor->coefficients[factor->coefficient_count - 1],
                                generator_minimal_polynomial), generator_minimal_polynomial);
                    ++factor_count;
                }
                break;
            }
            k = integer_subtract(local_stack, output_stack, k, &one);
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}