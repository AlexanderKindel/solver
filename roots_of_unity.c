#include "declarations.h"

struct NestedPolynomial*double_generator_number_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct NestedPolynomial*a, struct NestedPolynomial*b,
    struct RationalPolynomial*inner_annulling_polynomial,
    struct NestedPolynomial*outer_annulling_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*out =
        nested_polynomial_euclidean_remainder(output_stack, local_stack,
            number_field_polynomial_multiply(local_stack, output_stack, a, b,
                inner_annulling_polynomial), outer_annulling_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void quick_sort(struct Stack*stack_a, struct Stack*stack_b, struct Number**roots,
    struct RationalInterval**argument_estimates, size_t start, size_t end)
{
    if (start < end)
    {
        size_t pivot_index = start;
        for (size_t i = start + 1; i < end; ++i)
        {
            if (rational_compare(stack_a, stack_b, argument_estimates[pivot_index]->min,
                argument_estimates[i]->max) > 0)
            {
                POINTER_SWAP(argument_estimates[i], argument_estimates[pivot_index]);
                POINTER_SWAP(roots[i], roots[pivot_index]);
                pivot_index = i;
            }
        }
        if (pivot_index > 1)
        {
            quick_sort(stack_a, stack_b, roots, argument_estimates, start, pivot_index);
        }
        quick_sort(stack_a, stack_b, roots, argument_estimates, pivot_index + 1, end);
    }
}

void roots_of_unity_sort(struct Stack*stack_a, struct Stack*stack_b, struct Number**roots,
    struct Integer*degree)
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t degree_size_t = integer_to_size_t(degree);
    struct Rational*interval_size = rational_integer_divide(stack_a, stack_b, pi.min, degree);
    struct RationalInterval**argument_estimates =
        ARRAY_ALLOCATE(stack_a, degree_size_t - 1, struct FloatInterval*);
    for (size_t i = 1; i < degree_size_t; ++i)
    {
        argument_estimates[i - 1] = number_rational_argument_estimate(stack_a, stack_b, roots[i],
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
    struct Number**out = roots_of_unity + integer_to_size_t(integer_half(stack_a,
        integer_multiply(stack_a, stack_b, degree, degree_minus_one)));
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
        if (degree_size_t % prime_size_t == 0)
        {
            size_t quotient = degree_size_t / prime_size_t;
            struct Number**prime_degree_roots = get_roots_of_unity(stack_a, stack_b, prime);
            struct Number**quotient_degree_roots =
                get_roots_of_unity(stack_a, stack_b, integer_from_size_t(stack_a, quotient));
            for (size_t i = 0; i < prime_size_t; ++i)
            {
                for (size_t j = 0; j < quotient; ++j)
                {
                    out[i] = number_multiply(&permanent_stack, stack_b, prime_degree_roots[i],
                        number_root(stack_a, stack_b, quotient_degree_roots[j], prime));
                }
            }
            roots_of_unity_sort(stack_a, stack_b, out, degree);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        ++prime_index;
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
            if (integer_equals(integer_euclidean_remainder(stack_a, stack_b,
                integer_exponentiate(stack_a, stack_b, generator,
                    integer_from_size_t(stack_a, degree_minus_one_size_t /
                        degree_minus_one_factors[i])), degree), &one))
            {
                generator = integer_add(stack_a, generator, &one);
                goto generator_not_found;
            }
        }
        break;
    }
    struct RationalPolynomial*n_minus_first_root_annulling_polynomial =
        polynomial_allocate(stack_a, degree_size_t);
    n_minus_first_root_annulling_polynomial->coefficients[0] =
        &(struct Rational) { &INT(1, -), &one };
    struct NestedPolynomial*nth_root_annulling_polynomial =
        polynomial_allocate(stack_a, degree_size_t + 1);
    nth_root_annulling_polynomial->coefficients[0] = polynomial_allocate(stack_a, 1);
    nth_root_annulling_polynomial->coefficients[0]->coefficients[0] =
        n_minus_first_root_annulling_polynomial->coefficients[0];
    for (size_t i = 1; i < degree_minus_one_size_t; ++i)
    {
        n_minus_first_root_annulling_polynomial->coefficients[i] = &rational_zero;
        nth_root_annulling_polynomial->coefficients[i] = polynomial_zero;
    }
    n_minus_first_root_annulling_polynomial->coefficients[degree_minus_one_size_t] = &rational_one;
    nth_root_annulling_polynomial->coefficients[degree_minus_one_size_t] = polynomial_zero;
    nth_root_annulling_polynomial->coefficients[degree_size_t] = &rational_polynomial_one;

    //The resolvents array, as well as resolvent_multiples_in_terms_of_degree_minus_first_roots and
    //resolvent_product_values, exists for the purpose of calculating the entries of the
    //resolvent_values array. Conceptually, there are degree_minus_one resolvents, but the
    //resolvent_value entry corresponding to one of them is known a priori, so it is excluded from
    //the resolvents array, and thus also from resolvent_count. resolvent_values puts the entry
    //known a priori at index 0 (out of necessity), so resolvents[i] corresponds to
    //resolvent_values[i + 1].
    size_t resolvent_count = degree_minus_one_size_t - 1;
    struct NestedPolynomial**resolvents =
        ARRAY_ALLOCATE(stack_a, resolvent_count, struct NestedPolynomial*);
    resolvents[0] = polynomial_allocate(stack_a, 2);
    resolvents[0]->coefficients[0] = polynomial_zero;
    resolvents[0]->coefficients[1] = &rational_polynomial_one;
    for (size_t i = 1; i < resolvent_count; ++i)
    {
        resolvents[i] = resolvents[0];
    }
    struct Integer*generator_power = &one;
    struct Rational***resolvent_multiples_in_terms_of_degree_minus_first_roots =
        ARRAY_ALLOCATE(stack_a, resolvent_count, struct Rational**);
    for (size_t i = 0; i < resolvent_count; ++i)
    {
        generator_power = integer_euclidean_remainder(stack_a, stack_b,
            integer_multiply(stack_a, stack_b, generator_power, generator), degree);
        size_t degree_minus_first_root_exponent = 1;
        resolvent_multiples_in_terms_of_degree_minus_first_roots[i] =
            ARRAY_ALLOCATE(stack_a, degree_minus_one_size_t, struct Rational*);
        for (size_t j = 0; j < resolvent_count; ++j)
        {
            resolvent_multiples_in_terms_of_degree_minus_first_roots[i][j] = &rational_zero;
            size_t generator_power_size_t = integer_to_size_t(generator_power);
            struct NestedPolynomial*resolvent_term =
                polynomial_allocate(stack_a, generator_power_size_t + 1);
            for (size_t k = 0; k < generator_power_size_t; ++k)
            {
                resolvent_term->coefficients[k] = polynomial_zero;
            }
            resolvent_term->coefficients[generator_power_size_t] =
                polynomial_allocate(stack_a, degree_minus_first_root_exponent + 1);
            for (size_t k = 0; k < degree_minus_first_root_exponent; ++k)
            {
                resolvent_term->coefficients[generator_power_size_t]->coefficients[k] =
                    &rational_zero;
            }
            resolvent_term->coefficients[generator_power_size_t]->coefficients
                [degree_minus_first_root_exponent] = &rational_one;
            resolvents[j] = nested_polynomial_add(stack_a, stack_b, resolvents[j], resolvent_term);
            degree_minus_first_root_exponent =
                (degree_minus_first_root_exponent + i + 1) % degree_minus_one_size_t;
        }
        resolvent_multiples_in_terms_of_degree_minus_first_roots[i][resolvent_count] =
            &rational_zero;
    }
    struct NestedPolynomial*resolvent_power = &nested_polynomial_one;
    for (size_t i = resolvent_count; i-- > 0;)
    {
        resolvent_power = double_generator_number_multiply(stack_a, stack_b, resolvent_power,
            resolvents[0], n_minus_first_root_annulling_polynomial,
            nth_root_annulling_polynomial);
        struct NestedPolynomial*resolvent_product =
            double_generator_number_multiply(stack_a, stack_b, resolvent_power, resolvents[i],
                n_minus_first_root_annulling_polynomial, nth_root_annulling_polynomial);
        if (resolvent_product->coefficient_count > 1)
        {
            for (size_t j = 0; j < resolvent_product->coefficients[1]->coefficient_count; ++j)
            {
                resolvent_multiples_in_terms_of_degree_minus_first_roots[i][j] =
                    rational_subtract(stack_a, stack_b,
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[i][j],
                        resolvent_product->coefficients[1]->coefficients[j]);
            }
        }
        if (resolvent_product->coefficient_count)
        {
            for (size_t j = 0; j < resolvent_product->coefficients[0]->coefficient_count; ++j)
            {
                resolvent_multiples_in_terms_of_degree_minus_first_roots[i][j] =
                    rational_add(stack_a, stack_b,
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[i][j],
                        resolvent_product->coefficients[0]->coefficients[j]);
            }
        }
    }
    struct Number**degree_minus_first_roots =
        get_roots_of_unity(stack_a, stack_b, degree_minus_one);
    struct Number**resolvent_product_values =
        ARRAY_ALLOCATE(stack_a, resolvent_count, struct Number*);
    for (size_t i = 0; i < resolvent_count; ++i)
    {
        resolvent_product_values[i] = number_rational_initialize(stack_a, &rational_zero);
        for (size_t j = 0; j < degree_minus_one_size_t; ++j)
        {
            resolvent_product_values[i] =
                number_add(stack_a, stack_b, resolvent_product_values[i],
                    number_rational_multiply(stack_a, stack_b, degree_minus_first_roots[j],
                        resolvent_multiples_in_terms_of_degree_minus_first_roots[i][j]));
        }
    }
    struct Number**resolvent_values =
        ARRAY_ALLOCATE(stack_a, degree_minus_one_size_t, struct Number*);
    resolvent_values[0] =
        number_rational_initialize(stack_a, &(struct Rational){&INT(1, -), &one});
    resolvent_values[1] =
        number_root(stack_a, stack_b, resolvent_product_values[0], degree_minus_one);
    for (size_t i = 1; i < resolvent_count; ++i)
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
        out[i] = number_rational_initialize(stack_a, &rational_zero);
        size_t degree_minus_first_root_exponent = 0;
        for (size_t j = 0; j < degree_minus_one_size_t; ++j)
        {
            out[i] = number_add(stack_a, stack_b, number_multiply(stack_a, stack_b,
                degree_minus_first_roots[degree_minus_first_root_exponent], resolvent_values[j]),
                out[i]);
            degree_minus_first_root_exponent =
                (degree_minus_first_root_exponent + i) % degree_minus_one_size_t;
        }
        out[i] = number_eliminate_linear_dependencies(&permanent_stack, stack_b,
            number_rational_multiply(stack_a, stack_b, out[i], degree_minus_one_reciprocal));
    }
    roots_of_unity_sort(stack_a, stack_b, out, degree);
    stack_a->cursor = stack_a_savepoint;
    return out;
}