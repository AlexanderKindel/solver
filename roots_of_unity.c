#include "declarations.h"

void quick_sort(struct Stack*stack_a, struct Stack*stack_b, struct Number**roots,
    struct FloatInterval**argument_estimates, size_t start, size_t end)
{
    if (start < end)
    {
        size_t final_pivot_index = 0;
        for (size_t i = 0; i < end; ++i)
        {
            if (float_compare(stack_a, stack_b, argument_estimates[i]->max,
                argument_estimates[end]->min) < 0)
            {
                POINTER_SWAP(argument_estimates[i], argument_estimates[final_pivot_index]);
                POINTER_SWAP(roots[i], roots[final_pivot_index]);
                ++i;
            }
        }
        POINTER_SWAP(argument_estimates[final_pivot_index], argument_estimates[end]);
        POINTER_SWAP(roots[final_pivot_index], roots[end]);
        if (final_pivot_index > 1)
        {
            quick_sort(stack_a, stack_b, roots, argument_estimates, start, final_pivot_index - 1);
        }
        quick_sort(stack_a, stack_b, roots, argument_estimates, final_pivot_index + 1, end);
    }
}

void roots_of_unity_sort(struct Stack*stack_a, struct Stack*stack_b, struct Number**roots,
    struct Integer*degree)
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t degree_size_t = integer_to_size_t(degree);
    struct Rational*interval_size = rational_integer_divide(stack_a, stack_b, pi.min, degree);
    struct FloatInterval**argument_estimates =
        ARRAY_ALLOCATE(stack_a, degree_size_t - 1, struct FloatInterval*);
    for (size_t i = 1; i < degree_size_t; ++i)
    {
        number_float_argument_estimate(stack_a, stack_b, argument_estimates[i - 1], roots[i],
            interval_size);
    }
    quick_sort(stack_a, stack_b, &roots[1], argument_estimates, 0, degree_size_t - 1);
    stack_a->cursor = stack_a_savepoint;
}

struct Number**get_roots_of_unity(struct Stack*stack_a, struct Stack*stack_b, struct Integer*degree)
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t degree_size_t = integer_to_size_t(degree);
    size_t degree_minus_one_size_t = degree_size_t - 1;
    struct Integer*degree_minus_one = integer_from_size_t(stack_a, degree_minus_one_size_t);
    struct Number**out = roots_of_unity + integer_to_size_t(
        integer_half(stack_a, integer_multiply(stack_a, stack_b, degree, degree_minus_one)));
    if (out[0])
    {
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    size_t degree_square_root = integer_to_size_t(integer_square_root(stack_a, stack_b, degree));
    size_t prime_index = 0;
    while (true)
    {
        struct Integer*prime = get_prime(stack_a, stack_b, prime_index);
        size_t prime_size_t = integer_to_size_t(prime);
        if (prime_size_t > degree_square_root)
        {
            break;
        }
        if (degree_size_t % prime_size_t != 0)
        {
            size_t quotient = degree_size_t / prime_size_t;
            struct Number**prime_degree_roots = get_roots_of_unity(stack_a, stack_b, prime);
            struct Number**quotient_degree_roots =
                get_roots_of_unity(stack_a, stack_b, integer_from_size_t(stack_a, quotient));
            for (size_t i = 0; i < prime_size_t; ++i)
            {
                for (size_t j = 0; j < quotient; ++j)
                {
                    out[i] = number_copy(&permanent_stack, number_multiply(stack_a, stack_b,
                        prime_degree_roots[i], number_exponentiate(stack_a, stack_b,
                            quotient_degree_roots[j], &(struct Rational){&one, prime})));
                }
            }
            roots_of_unity_sort(stack_a, stack_b, out, degree);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
    }
    size_t*degree_minus_one_factors;
    size_t factor_count =
        size_t_factor(stack_a, stack_b, &degree_minus_one_factors, degree_minus_one_size_t);
    struct Integer*generator = &INT(2, +);
    while (true)
    {
    generator_not_found:
        for (size_t i = 0; i < factor_count; ++i)
        {
            if (!integer_equals(integer_euclidean_remainder(stack_a, stack_b,
                integer_exponentiate(stack_a, stack_b, generator, integer_from_size_t(stack_a,
                    degree_minus_one_size_t / degree_minus_one_factors[i])), degree), &one))
            {
                generator = integer_add(stack_a, generator, &one);
                goto generator_not_found;
            }
        }
        break;
    }
    struct RationalPolynomial*resolvent_generator_annulling_polynomials[2] =
        { polynomial_allocate(stack_a, degree_size_t),
            polynomial_allocate(stack_a, degree_size_t + 1) };
    resolvent_generator_annulling_polynomials[0]->coefficients[0] =
        &(struct Rational) { &INT(1, -), &one };
    resolvent_generator_annulling_polynomials[1]->coefficients[0] =
        resolvent_generator_annulling_polynomials[0]->coefficients[0];
    for (size_t i = 1; i < degree_minus_one_size_t; ++i)
    {
        resolvent_generator_annulling_polynomials[0]->coefficients[i] = &rational_zero;
        resolvent_generator_annulling_polynomials[1]->coefficients[i] = &rational_zero;
    }
    resolvent_generator_annulling_polynomials[0]->coefficients[degree_minus_one_size_t] =
        &rational_one;
    resolvent_generator_annulling_polynomials[1]->coefficients[degree_minus_one_size_t] =
        &rational_zero;
    resolvent_generator_annulling_polynomials[1]->coefficients[degree_size_t] = &rational_one;
    struct AlgebraicNumber**resolvents =
        ARRAY_ALLOCATE(stack_a, degree_minus_one_size_t - 1, struct AlgebraicNumber*);
    size_t resolvent_count_minus_one = degree_minus_one_size_t - 1;
    for (size_t i = 0; i < resolvent_count_minus_one; ++i)
    {
        resolvents[i] = ALLOCATE(stack_a, struct AlgebraicNumber);
        resolvents[i]->next_term = 0;
        resolvents[i]->term_coefficient = &rational_one;
        resolvents[i]->generator_degrees[0] = 0;
        resolvents[i]->generator_degrees[1] = 1;
    }
    struct Integer*generator_power = &one;
    struct Rational**resolvent_multiples_in_terms_of_degree_minus_first_roots =
        ARRAY_ALLOCATE(stack_a, resolvent_count_minus_one * degree_minus_one_size_t,
            struct Rational*);
    for (size_t i = 0; i < resolvent_count_minus_one; ++i)
    {
        size_t resolvent_multiple_index = degree_minus_one_size_t * i;
        resolvent_multiples_in_terms_of_degree_minus_first_roots[resolvent_multiple_index] =
            &rational_zero;
        generator_power = integer_euclidean_remainder(stack_a, stack_b,
            integer_multiply(stack_a, stack_b, generator_power, generator), degree);
        size_t degree_minus_first_root_exponent = 1;
        for (size_t j = 1; j < degree_minus_one_size_t; ++j)
        {
            struct AlgebraicNumber*resolvent_term = ALLOCATE(stack_a, struct AlgebraicNumber);
            resolvent_term->next_term = 0;
            resolvent_term->term_coefficient = &rational_one;
            resolvent_term->generator_degrees[0] = degree_minus_first_root_exponent;
            resolvent_term->generator_degrees[1] = integer_to_size_t(generator_power);
            resolvents[j - 1] =
                algebraic_number_add(stack_a, stack_b, resolvents[j - 1], resolvent_term);
            resolvent_multiples_in_terms_of_degree_minus_first_roots[resolvent_multiple_index + j] =
                &rational_zero;
            degree_minus_first_root_exponent =
                (degree_minus_first_root_exponent + i + 1) % degree_minus_one_size_t;
        }
    }
    struct AlgebraicNumber*resolvent_power = ALLOCATE(stack_a, struct AlgebraicNumber);
    resolvent_power->next_term = 0;
    resolvent_power->term_coefficient = &rational_one;
    resolvent_power->generator_degrees[0] = 0;
    resolvent_power->generator_degrees[1] = 0;
    for (size_t i = resolvent_count_minus_one; i-- > 0;)
    {
        resolvent_power = algebraic_number_multiply(stack_a, stack_b, resolvent_power,
            resolvents[0], resolvent_generator_annulling_polynomials);
        struct AlgebraicNumber*resolvent_product = algebraic_number_multiply(stack_a, stack_b,
            resolvent_power, resolvents[i], resolvent_generator_annulling_polynomials);
        struct AlgebraicNumber*resolvent_product_term = resolvent_product;
        while (resolvent_product_term)
        {
            switch (resolvent_product_term->generator_degrees[1])
            {
            case 0:
            {
                size_t index =
                    degree_minus_one_size_t * i + resolvent_product_term->generator_degrees[0];
                resolvent_multiples_in_terms_of_degree_minus_first_roots[index] =
                    rational_add(stack_a, stack_b, resolvent_product_term->term_coefficient,
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[index]);
            }
            case 1:
            {
                size_t index =
                    degree_minus_one_size_t * i + resolvent_product_term->generator_degrees[0];
                resolvent_multiples_in_terms_of_degree_minus_first_roots[index] =
                    rational_subtract(stack_a, stack_b, resolvent_product_term->term_coefficient,
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[index]);
            }
            }
            resolvent_product_term = resolvent_product_term->next_term;
        }
    }
    struct Number**degree_minus_first_roots =
        get_roots_of_unity(stack_a, stack_b, degree_minus_one);
    struct Number**resolvent_product_values =
        ARRAY_ALLOCATE(stack_a, resolvent_count_minus_one, struct Number*);
    for (size_t i = 0; i < resolvent_count_minus_one; ++i)
    {
        size_t resolvent_multiple_index = degree_minus_one_size_t * i;
        resolvent_product_values[i] = roots_of_unity[0];
        for (size_t j = 0; j < degree_minus_one_size_t; ++j)
        {
            resolvent_product_values[i] =
                number_add(stack_a, stack_b, resolvent_product_values[i],
                    number_rational_multiply(stack_a, stack_b, degree_minus_first_roots[i],
                        resolvent_multiples_in_terms_of_degree_minus_first_roots
                            [resolvent_multiple_index]));
        }
    }
    struct Number**resolvent_values =
        ARRAY_ALLOCATE(stack_a, degree_minus_one_size_t, struct Number*);
    resolvent_values[0] = ALLOCATE(stack_a, struct Number);
    resolvent_values[0]->operation = 'r';
    resolvent_values[0]->value.numerator = integer_initialize(stack_a, 1, -1);
    resolvent_values[0]->value.denominator = &one;
    resolvent_values[1] = number_exponentiate(stack_a, stack_b,
        resolvent_product_values[0], &(struct Rational){&one, degree_minus_one});
    for (size_t i = 1; i < resolvent_count_minus_one; ++i)
    {
        resolvent_values[i + 1] = number_divide(stack_a, stack_b, resolvent_values[1],
            resolvent_product_values[i]);
    }
    struct Rational*degree_minus_one_reciprocal = ALLOCATE(stack_a, struct Rational);
    degree_minus_one_reciprocal->numerator = &one;
    degree_minus_one_reciprocal->denominator = degree_minus_one;
    out[0] = roots_of_unity[0];
    for (size_t i = 1; i < degree_size_t; ++i)
    {
        out[i] = ALLOCATE(stack_a, struct Number);
        out[i]->operation = 'r';
        out[i]->value = rational_zero;
        size_t degree_minus_first_root_exponent = 0;
        for (size_t j = 0; j < degree_minus_one_size_t; ++j)
        {
            out[i] = number_add(stack_a, stack_b, number_multiply(stack_a, stack_b,
                degree_minus_first_roots[degree_minus_first_root_exponent], resolvent_values[j]),
                out[i]);
            degree_minus_first_root_exponent =
                (degree_minus_first_root_exponent + i) % degree_minus_one_size_t;
        }
        out[i] = number_copy(&permanent_stack,
            number_rational_multiply(stack_a, stack_b, out[i], degree_minus_one_reciprocal));
    }
    roots_of_unity_sort(stack_a, stack_b, out, degree);
    stack_a->cursor = stack_a_savepoint;
    return out;
}