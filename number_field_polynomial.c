#include "declarations.h"

struct NestedPolynomial*number_field_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    return polynomial_multiply(&number_field_element_operations.ring_operations, output_stack,
        local_stack, (struct Polynomial*)a, (struct Polynomial*)b, generator_minimal_polynomial);
}

struct NestedPolynomial*number_field_polynomial_element_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct RationalPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    return polynomial_multiply_by_coefficient(&number_field_element_operations.ring_operations,
        output_stack, local_stack, (struct Polynomial*)a, b, generator_minimal_polynomial);
}

void number_field_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct NestedPolynomial*dividend,
    struct NestedPolynomial*divisor, struct RationalPolynomial*generator_minimal_polynomial)
{
    field_polynomial_euclidean_divide(&number_field_element_operations, output_stack, local_stack,
        out, (struct Polynomial*)dividend, (struct Polynomial*)divisor,
        generator_minimal_polynomial);
}

struct NestedPolynomial*number_field_polynomial_gcd(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    return generic_gcd(&number_field_polynomial_operations, output_stack, local_stack, a, b,
        generator_minimal_polynomial);
}

size_t number_field_polynomial_squarefree_factor(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial**out,
    struct RationalPolynomial*generator_minimal_polynomial)
{
    return polynomial_squarefree_factor(&number_field_polynomial_operations,
        nested_polynomial_derivative, output_stack, local_stack, (struct Polynomial*)a,
        (struct Polynomial**)out, generator_minimal_polynomial);
}

size_t number_field_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct NestedPolynomial*a, struct RationalPolynomial*generator_minimal_polynomial,
    struct NestedPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial**squarefree_factors = stack_slot_allocate(local_stack,
        (a->coefficient_count - 1) * sizeof(struct NestedPolynomial*),
        _Alignof(struct NestedPolynomial*));
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
            struct NestedPolynomial*e = (struct NestedPolynomial*)&polynomial_zero;
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
            if (rational_polynomial_gcd(local_stack, output_stack,
                rational_polynomial_derivative(local_stack, output_stack, resultant),
                resultant)->coefficient_count < 2)
            {
                struct IntegerPolynomial**resultant_factors = stack_slot_allocate(local_stack,
                    (resultant->coefficient_count - 1) * sizeof(struct IntegerPolynomial*),
                    _Alignof(struct IntegerPolynomial*));
                size_t resultant_factor_count =
                    primitive_integer_polynomial_factor(local_stack, output_stack,
                        rational_polynomial_primitive_part(local_stack, output_stack, resultant),
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
                    e = (struct NestedPolynomial*)&polynomial_zero;
                    for (size_t k = 0; k < resultant_factors[j]->coefficient_count; ++k)
                    {
                        struct RationalPolynomial*coefficient = polynomial_allocate(local_stack, 1);
                        coefficient->coefficients[0] =
                            &(struct Rational) { resultant_factors[j]->coefficients[k], &one };
                        e = nested_polynomial_add(local_stack, output_stack, e,
                            nested_polynomial_rational_polynomial_multiply(local_stack,
                                output_stack, power, coefficient));
                        power = number_field_polynomial_multiply(local_stack, output_stack, power,
                            d, generator_minimal_polynomial);
                    }
                    struct NestedPolynomial*factor = number_field_polynomial_gcd(local_stack,
                        output_stack, e, squarefree_factors[i], generator_minimal_polynomial);
                    out[factor_count] =
                        number_field_polynomial_element_multiply(output_stack, local_stack, factor,
                            number_field_element_reciprocal(local_stack, output_stack,
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