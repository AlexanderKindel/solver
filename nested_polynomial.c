#include "declarations.h"

struct NestedPolynomial*nested_polynomial_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor)
{
    struct NestedPolynomialDivision division;
    nested_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return division.remainder;
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
    struct RationalPolynomial*a_content = nested_polynomial_content(local_stack, output_stack, a);
    struct RationalPolynomial*b_content = nested_polynomial_content(local_stack, output_stack, b);
    struct RationalPolynomial*t = rational_polynomial_multiply(local_stack, output_stack,
        rational_polynomial_exponentiate(local_stack, output_stack, a_content,
            size_t_to_integer(local_stack, b->coefficient_count - 1)),
        rational_polynomial_exponentiate(local_stack, output_stack, b_content,
            size_t_to_integer(local_stack, a->coefficient_count - 1)));
    a = nested_polynomial_rational_polynomial_divide(local_stack, output_stack, a, a_content);
    b = nested_polynomial_rational_polynomial_divide(local_stack, output_stack, b, b_content);
    struct RationalPolynomial*g = &rational_polynomial_one;
    struct RationalPolynomial*h = &rational_polynomial_one;
    if (b->coefficient_count > a->coefficient_count)
    {
        POINTER_SWAP(a, b);
    }
    while (b->coefficient_count > 1)
    {
        size_t degree_size_t = a->coefficient_count - b->coefficient_count;
        struct Integer*degree = size_t_to_integer(local_stack, degree_size_t);
        struct NestedPolynomial*remainder =
            nested_polynomial_get_remainder(local_stack, output_stack,
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
            h = rational_polynomial_get_quotient(local_stack, output_stack, term,
                rational_polynomial_exponentiate(local_stack, output_stack, h,
                    size_t_to_integer(local_stack, degree_size_t - 1)));
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
        out = rational_polynomial_get_quotient(local_stack, output_stack, term,
            rational_polynomial_exponentiate(local_stack, output_stack, h,
                size_t_to_integer(local_stack, a->coefficient_count - 2)));
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

struct RationalPolynomial*number_field_element_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = rational_polynomial_get_remainder(output_stack, local_stack,
        rational_polynomial_multiply(local_stack, output_stack, a, b),
        generator_minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct PolynomialExtendedGCDInfo info;
    rational_polynomial_get_extended_gcd(local_stack, output_stack, &info, a,
        generator_minimal_polynomial);
    struct RationalPolynomial*out = rational_polynomial_rational_multiply(output_stack, local_stack,
        (struct RationalPolynomial*)info.a_coefficient,
        rational_get_reciprocal(local_stack,
            ((struct RationalPolynomial*)info.gcd)->coefficients[0]));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_field_element_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct RationalPolynomial*dividend,
    struct RationalPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*out = number_field_element_multiply(output_stack, local_stack,
        dividend, number_field_element_get_reciprocal(local_stack, output_stack, divisor,
            generator_minimal_polynomial), generator_minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t number_field_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct NestedPolynomial*a,
    struct RationalPolynomial*generator_minimal_polynomial, struct NestedPolynomial**out)
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
            struct NestedPolynomial*power = &nested_polynomial_one;
            struct NestedPolynomial*d = polynomial_allocate(local_stack, 2);
            d->coefficients[0] = polynomial_allocate(local_stack, 2);
            d->coefficients[0]->coefficients[0] = &rational_zero;
            d->coefficients[0]->coefficients[1] = &rational_one;
            d->coefficients[1] = polynomial_allocate(local_stack, 1);
            d->coefficients[1]->coefficients[0] = &(struct Rational) { k, &one };
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
                rational_polynomial_derivative(local_stack, output_stack, resultant),
                resultant)->coefficient_count < 2)
            {
                struct IntegerPolynomial**resultant_factors = ARRAY_ALLOCATE(local_stack,
                    resultant->coefficient_count - 1, struct IntegerPolynomial*);
                size_t resultant_factor_count =
                    primitive_integer_polynomial_factor(local_stack, output_stack,
                        rational_polynomial_get_primitive_part(local_stack, output_stack,
                            resultant),
                        resultant_factors);
                k->sign *= -1;
                d = polynomial_allocate(local_stack, 2);
                d->coefficients[0] = polynomial_allocate(local_stack, 2);
                d->coefficients[0]->coefficients[0] = &rational_zero;
                d->coefficients[0]->coefficients[1] = &(struct Rational) { k, &one };
                d->coefficients[1] = &rational_polynomial_one;
                for (size_t j = 0; j < resultant_factor_count; ++j)
                {
                    power = &nested_polynomial_one;
                    e = polynomial_zero;
                    for (size_t k = 0; k < resultant_factors[j]->coefficient_count; ++k)
                    {
                        e = nested_polynomial_add(local_stack, output_stack, e,
                            nested_polynomial_rational_polynomial_multiply(local_stack,
                                output_stack, power,
                                POLY(Rational, 1,
                                    &(struct Rational) {resultant_factors[j]->coefficients[k], &one })));
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