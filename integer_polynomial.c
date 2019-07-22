#include "declarations.h"

struct IntegerPolynomial*integer_polynomial_copy(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    return polynomial_copy(integer_copy_to_stack, output_stack, (struct Polynomial*)a);
}

bool integer_polynomial_equals(struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    return polynomial_equals(&integer_operations.ring_operations, (struct Polynomial*)a,
        (struct Polynomial*)b);
}

struct IntegerPolynomial*integer_polynomial_add(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    return polynomial_add(&integer_operations.ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b);
}

struct IntegerPolynomial*integer_polynomial_negative(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    return polynomial_negative(&integer_operations.ring_operations, output_stack,
        (struct Polynomial*)a);
}

struct IntegerPolynomial*integer_polynomial_subtract(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*minuend, struct IntegerPolynomial*subtrahend)
{
    return polynomial_subtract(&integer_operations.ring_operations, output_stack, local_stack,
        (struct Polynomial*)minuend, (struct Polynomial*)subtrahend);
}

struct IntegerPolynomial*integer_polynomial_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    return polynomial_multiply(&integer_operations.ring_operations, output_stack, local_stack,
        (struct Polynomial*)a, (struct Polynomial*)b, 0);
}

struct IntegerPolynomial*integer_polynomial_integer_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Integer*b)
{
    return polynomial_multiply_by_coefficient(&integer_operations.ring_operations, output_stack,
        local_stack, (struct Polynomial*)a, b, 0);
}

void integer_polynomial_euclidean_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialDivision*out, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor)
{
    polynomial_euclidean_divide(&integer_operations, output_stack, local_stack, out,
        (struct Polynomial*)dividend, (struct Polynomial*)divisor, 0);
}

struct IntegerPolynomial*integer_polynomial_euclidean_quotient(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor)
{
    struct PolynomialDivision division;
    integer_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct IntegerPolynomial*)division.quotient;
}

struct IntegerPolynomial*integer_polynomial_euclidean_remainder(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor)
{
    struct PolynomialDivision division;
    integer_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return (struct IntegerPolynomial*)division.remainder;
}

struct IntegerPolynomial*integer_polynomial_integer_divide(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*dividend, struct Integer*divisor)
{
    return polynomial_divide_by_coefficient(integer_euclidean_quotient, output_stack, local_stack,
        (struct Polynomial*)dividend, divisor);
}

struct Integer*integer_polynomial_content(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a)
{
    return polynomial_content(&integer_operations, output_stack, local_stack,
        (struct Polynomial*)a, 0);
}

struct IntegerPolynomial*integer_polynomial_primitive_part(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_integer_divide(output_stack, local_stack, a,
        integer_polynomial_content(local_stack, output_stack, a));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    struct PolynomialExtendedGCDInfo info;
    integer_polynomial_extended_gcd(output_stack, local_stack, &info, a, b);
    return (struct IntegerPolynomial*)info.gcd;
}

void integer_polynomial_extended_gcd(struct Stack*output_stack, struct Stack*local_stack,
    struct PolynomialExtendedGCDInfo*out, struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    if (b->coefficient_count > a->coefficient_count)
    {
        POINTER_SWAP(a, b);
    }
    if (b->coefficient_count == 0)
    {
        out->gcd = (struct Polynomial*)integer_polynomial_copy(output_stack, a);
        out->a_coefficient = (struct Polynomial*)&integer_polynomial_one;
        return;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*a_content = integer_polynomial_content(local_stack, output_stack, a);
    struct Integer*b_content = integer_polynomial_content(local_stack, output_stack, b);
    struct Integer*d = integer_gcd(local_stack, output_stack, a_content, b_content);
    struct Integer*g = &one;
    struct Integer*h = &one;
    struct IntegerPolynomial*previous_a_coefficient = &integer_polynomial_one;
    struct IntegerPolynomial*previous_gcd_multiple =
        integer_polynomial_integer_divide(local_stack, output_stack, a, a_content);
    out->a_coefficient = &polynomial_zero;
    out->gcd = (struct Polynomial*)integer_polynomial_integer_divide(local_stack, output_stack, b,
        b_content);
    while (true)
    {
        struct Integer*degree = integer_from_size_t(local_stack,
            previous_gcd_multiple->coefficient_count - out->gcd->coefficient_count);
        struct Integer*multiple = integer_exponentiate(local_stack, output_stack,
            out->gcd->coefficients[out->gcd->coefficient_count - 1],
            integer_add(local_stack, degree, &one));
        struct PolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, previous_gcd_multiple,
                multiple), (struct IntegerPolynomial*)out->gcd);
        switch (division.remainder->coefficient_count)
        {
        case 0:
        {
            struct Integer*gcd_content = integer_polynomial_content(local_stack, output_stack,
                (struct IntegerPolynomial*)out->gcd);
            out->gcd = (struct Polynomial*)integer_polynomial_integer_multiply(output_stack,
                local_stack, integer_polynomial_integer_divide(local_stack, output_stack,
                    (struct IntegerPolynomial*)out->gcd, gcd_content), d);
            out->a_coefficient =
                (struct Polynomial*)integer_polynomial_integer_multiply(output_stack, local_stack,
                    integer_polynomial_integer_divide(local_stack, output_stack,
                        (struct IntegerPolynomial*)out->a_coefficient, gcd_content), d);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        case 1:
        {
            out->gcd = stack_polynomial_allocate(output_stack, 1);
            out->gcd->coefficients[0] = integer_copy_to_stack(output_stack, d);
            out->a_coefficient =
                (struct Polynomial*)integer_polynomial_integer_multiply(output_stack, local_stack,
                    (struct IntegerPolynomial*)out->a_coefficient, d);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        }
        struct Integer*divisor = integer_multiply(local_stack, output_stack, g,
            integer_exponentiate(local_stack, output_stack, h, degree));
        struct IntegerPolynomial*next_a_coefficient =
            integer_polynomial_subtract(local_stack, output_stack,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    previous_a_coefficient, multiple),
                integer_polynomial_multiply(local_stack, output_stack,
                    (struct IntegerPolynomial*)division.quotient,
                    (struct IntegerPolynomial*)out->a_coefficient));
        previous_a_coefficient = (struct IntegerPolynomial*)out->a_coefficient;
        out->a_coefficient = (struct Polynomial*)integer_polynomial_integer_divide(local_stack,
            output_stack, next_a_coefficient, divisor);
        previous_gcd_multiple = (struct IntegerPolynomial*)out->gcd;
        out->gcd = (struct Polynomial*)integer_polynomial_integer_divide(local_stack, output_stack,
            (struct IntegerPolynomial*)division.remainder, divisor);
        g = previous_gcd_multiple->coefficients[previous_gcd_multiple->coefficient_count - 1];
        struct Integer*h_exponent = integer_subtract(local_stack, output_stack, &one, degree);
        if (h_exponent->sign >= 0)
        {
            h = integer_multiply(local_stack, output_stack,
                integer_exponentiate(local_stack, output_stack, h, h_exponent),
                integer_exponentiate(local_stack, output_stack, g, degree));
        }
        else
        {
            h_exponent->sign = -h_exponent->sign;
            h = integer_euclidean_quotient(local_stack, output_stack,
                integer_exponentiate(local_stack, output_stack, g, degree),
                integer_exponentiate(local_stack, output_stack, h, h_exponent));
        }
    }
}

struct IntegerPolynomial*bound_coefficients(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct Integer*characteristic_power)
{
    struct IntegerPolynomial*out = stack_polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        struct Integer*remainder = integer_euclidean_remainder(local_stack, output_stack,
            a->coefficients[i], characteristic_power);
        if (remainder->sign < 0)
        {
            remainder = integer_add(local_stack, remainder, characteristic_power);
        }
        if (integer_compare(output_stack, local_stack, integer_doubled(local_stack, remainder),
            characteristic_power) > 0)
        {
            a->coefficients[i] =
                integer_subtract(output_stack, local_stack, remainder, characteristic_power);
        }
        else
        {
            a->coefficients[i] = integer_copy_to_stack(output_stack, remainder);
        }
    }
    return out;
}

struct IntegerPolynomial*find_valid_combination(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*combination, size_t combination_size,
    struct IntegerPolynomial**factor_candidates_start,
    struct IntegerPolynomial**factor_candidates_end, struct IntegerPolynomial**a,
    struct Integer*characteristic_power)
{
    if (combination_size == 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct IntegerPolynomial*bounded_combination =
            bound_coefficients(local_stack, output_stack, combination, characteristic_power);
        struct PolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, *a,
            (*a)->coefficients[(*a)->coefficient_count - 1]), bounded_combination);
        if (division.remainder && division.remainder->coefficient_count == 0)
        {
            do
            {
                *a = (struct IntegerPolynomial*)division.quotient;
                integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
                    integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                    (*a)->coefficients[(*a)->coefficient_count - 1]),
                    bounded_combination);
            } while (division.remainder && division.remainder->coefficient_count == 0);
            struct IntegerPolynomial*factor =
                integer_polynomial_primitive_part(output_stack, local_stack, bounded_combination);
            local_stack->cursor = local_stack_savepoint;
            return factor;
        }
        local_stack->cursor = local_stack_savepoint;
    }
    else
    {
        for (struct IntegerPolynomial**candidate = factor_candidates_start;
            *candidate <= *factor_candidates_end - combination_size; ++*candidate)
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
                integer_polynomial_multiply(local_stack, output_stack, combination, *candidate),
                combination_size - 1, factor_candidates_start + 1, factor_candidates_end, a,
                characteristic_power);
            local_stack->cursor = local_stack_savepoint;
            if (factor)
            {
                *factor_candidates_end -= 1;
                *candidate = *factor_candidates_end;
                return factor;
            }
        }
    }
    return 0;
}

size_t squarefree_integer_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*modded_a;
    struct IntegerPolynomial*gcd;
    size_t prime_index = 0;
    struct Integer*prime = primes[0];
    while (true)
    {
        modded_a = modded_polynomial_reduced(local_stack, output_stack, a, prime);
        if (modded_a->coefficient_count == a->coefficient_count)
        {
            gcd = modded_polynomial_gcd(local_stack, output_stack, modded_a,
                modded_polynomial_reduced(local_stack, output_stack,
                    integer_polynomial_derivative(local_stack, output_stack, modded_a), prime),
                prime);
            if (gcd->coefficient_count == 1)
            {
                break;
            }
        }
        ++prime_index;
        prime = get_prime(output_stack, local_stack, prime_index);
    }
    modded_a = modded_polynomial_monic(local_stack, output_stack, modded_a, prime);
    struct IntegerPolynomial**modded_a_factors = stack_slot_allocate(local_stack,
        (modded_a->coefficient_count - 1) * sizeof(struct IntegerPolynomial*),
        _Alignof(struct IntegerPolynomial*));
    size_t modded_a_factor_count = squarefree_modded_polynomial_factor(local_stack, output_stack,
        modded_a, prime, modded_a_factors);
    struct Integer*coefficient_bound = &zero;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        coefficient_bound = integer_add(local_stack, coefficient_bound,
            integer_multiply(local_stack, output_stack, a->coefficients[i], a->coefficients[i]));
    }
    struct Integer*a_degree_minus_one = integer_from_size_t(local_stack, a->coefficient_count - 2);
    struct Integer*k =
        integer_euclidean_quotient(local_stack, output_stack, a_degree_minus_one, &INT(2, +));
    struct Integer*leading_coefficient_magnitude =
        integer_magnitude(local_stack, a->coefficients[a->coefficient_count - 1]);
    coefficient_bound = integer_doubled(local_stack, integer_multiply(local_stack, output_stack,
        leading_coefficient_magnitude, integer_add(local_stack,
            integer_multiply(local_stack, output_stack,
                integer_square_root(local_stack, output_stack, coefficient_bound),
                n_choose_k(local_stack, output_stack, a_degree_minus_one, k)),
            integer_multiply(local_stack, output_stack, leading_coefficient_magnitude,
                n_choose_k(local_stack, output_stack, a_degree_minus_one,
                    integer_add(local_stack, k, &INT(1, -)))))));
    struct Integer*e = &one;
    struct Integer*prime_power = prime;
    while (integer_compare(output_stack, local_stack, prime_power, coefficient_bound) < 0)
    {
        e = integer_add(local_stack, e, &one);
        prime_power = integer_multiply(local_stack, output_stack, prime_power, prime);
    }
    struct IntegerPolynomial*product_of_unlifted_factors = modded_a;
    struct IntegerPolynomial*unmodded_b_times_c = modded_polynomial_monic(local_stack, output_stack,
        modded_polynomial_reduced(local_stack, output_stack, a, prime_power), prime_power);
    struct IntegerPolynomial**lifted_factors = stack_slot_allocate(local_stack,
        modded_a_factor_count * sizeof(struct IntegerPolynomial*),
        _Alignof(struct IntegerPolynomial*));
    size_t lifted_factor_count = 0;
    modded_a_factor_count -= 1;
    for (size_t i = 0; i < modded_a_factor_count; ++i)
    {
        product_of_unlifted_factors = modded_polynomial_euclidean_quotient(local_stack,
            output_stack, product_of_unlifted_factors, modded_a_factors[i], prime);
        struct Integer*lift_characteristic = prime;
        struct IntegerPolynomial*b = modded_a_factors[i];
        struct IntegerPolynomial*c = product_of_unlifted_factors;
        while (integer_compare(output_stack, local_stack, lift_characteristic, prime_power) < 0)
        {
            struct IntegerPolynomial*b_mod_prime =
                modded_polynomial_reduced(local_stack, output_stack, b, prime);
            struct IntegerPolynomial*c_mod_prime =
                modded_polynomial_reduced(local_stack, output_stack, c, prime);
            struct ExtendedGCDInfo gcd_info;
            modded_polynomial_extended_gcd(local_stack, output_stack, &gcd_info, b_mod_prime,
                c_mod_prime, prime);
            struct Integer*reciprocal = modded_integer_reciprocal(local_stack, output_stack,
                ((struct IntegerPolynomial*)gcd_info.gcd)->coefficients[0], prime);
            gcd_info.a_coefficient = modded_polynomial_multiply_by_coefficient(local_stack,
                output_stack, (struct IntegerPolynomial*)gcd_info.a_coefficient, reciprocal, prime);
            gcd_info.b_coefficient = modded_polynomial_multiply_by_coefficient(local_stack,
                output_stack, (struct IntegerPolynomial*)gcd_info.b_coefficient, reciprocal, prime);
            struct IntegerPolynomial*f = modded_polynomial_reduced(local_stack, output_stack,
                integer_polynomial_integer_divide(local_stack, output_stack,
                    integer_polynomial_subtract(local_stack, output_stack, unmodded_b_times_c,
                        integer_polynomial_multiply(local_stack, output_stack, b, c)),
                    lift_characteristic), prime);
            struct PolynomialDivision division;
            modded_polynomial_euclidean_divide(local_stack, output_stack, &division,
                modded_polynomial_multiply(local_stack, output_stack, f,
                (struct IntegerPolynomial*)gcd_info.b_coefficient, prime), b_mod_prime, prime);
            b = integer_polynomial_add(local_stack, output_stack, b,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                (struct IntegerPolynomial*)division.remainder, lift_characteristic));
            c = integer_polynomial_add(local_stack, output_stack, c,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    modded_polynomial_add(local_stack, output_stack,
                        modded_polynomial_multiply(local_stack, output_stack,
                        (struct IntegerPolynomial*)gcd_info.a_coefficient, f, prime),
                        modded_polynomial_multiply(local_stack, output_stack, c_mod_prime,
                        (struct IntegerPolynomial*)division.quotient, prime), prime),
                    lift_characteristic));
            lift_characteristic =
                integer_multiply(local_stack, output_stack, lift_characteristic, prime);
        }
        lifted_factors[lifted_factor_count] = b;
        ++lifted_factor_count;
        unmodded_b_times_c = c;
    }
    lifted_factors[lifted_factor_count] = unmodded_b_times_c;
    ++lifted_factor_count;
    struct IntegerPolynomial*a_times_leading_coefficient =
        integer_polynomial_integer_multiply(local_stack, output_stack, a,
            a->coefficients[a->coefficient_count - 1]);
    size_t factor_count = 0;
    size_t combination_size = 1;
    struct IntegerPolynomial*unfactored_component_of_a = a;
    struct IntegerPolynomial**lifted_factors_end = lifted_factors + lifted_factor_count;
    while (2 * combination_size < lifted_factor_count)
    {
        while (true)
        {
            struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
                &integer_polynomial_one, combination_size, lifted_factors, lifted_factors_end,
                &unfactored_component_of_a, prime_power);
            if (!factor)
            {
                break;
            }
            out[factor_count] = factor;
            ++factor_count;
        }
        ++combination_size;
    }
    if (2 * combination_size == lifted_factor_count)
    {
        lifted_factors_end -= 1;
        struct IntegerPolynomial*factor = find_valid_combination(output_stack, local_stack,
            *lifted_factors_end, combination_size - 1, lifted_factors, lifted_factors_end,
            &unfactored_component_of_a, prime_power);
        if (factor)
        {
            out[factor_count] = factor;
            ++factor_count;
        }
    }
    struct IntegerPolynomial*final_factor = stack_polynomial_allocate(local_stack, 1);
    final_factor->coefficients[0] =
        unfactored_component_of_a->coefficients[unfactored_component_of_a->coefficient_count - 1];
    for (struct IntegerPolynomial**factor = lifted_factors; factor < lifted_factors_end;
        factor += 1)
    {
        final_factor =
            integer_polynomial_multiply(local_stack, output_stack, final_factor, *factor);
    }
    if (final_factor->coefficient_count > 1)
    {
        out[factor_count] = integer_polynomial_primitive_part(output_stack, local_stack,
            bound_coefficients(local_stack, output_stack, final_factor, prime_power));
        ++factor_count;
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

bool list_contains_a(struct IntegerPolynomial**list, size_t list_length,
    struct IntegerPolynomial*a)
{
    for (size_t i = 0; i < list_length; ++i)
    {
        if (a->coefficient_count == list[i]->coefficient_count)
        {
            for (size_t j = 0; j < a->coefficient_count; ++j)
            {
                if (!integer_equals(a->coefficients[j], list[i]->coefficients[j]))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

size_t integer_polynomial_squarefree_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    return polynomial_squarefree_factor(&integer_polynomial_operations,
        integer_polynomial_derivative, output_stack, local_stack, (struct Polynomial*)a,
        (struct Polynomial**)out, 0);
}

size_t primitive_integer_polynomial_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    if (!a->coefficients[0]->value_count)
    {
        out[0] = stack_polynomial_allocate(output_stack, 2);
        out[0]->coefficients[0] = &zero;
        out[0]->coefficients[1] = &one;
        factor_count = 1;
        do
        {
            a = integer_polynomial_euclidean_quotient(local_stack, output_stack, a, *out);
        } while (!a->coefficients[0]->value_count);
    }
    if (a->coefficients[0]->value_count < 2)
    {
        local_stack->cursor = local_stack_savepoint;
        return factor_count;
    }
    struct IntegerPolynomial**squarefree_factors = stack_slot_allocate(local_stack,
        (a->coefficient_count - 1) * sizeof(struct IntegerPolynomial*),
        _Alignof(struct IntegerPolynomial*));
    size_t squarefree_factor_count =
        integer_polynomial_squarefree_factor(local_stack, output_stack, a, squarefree_factors);
    for (size_t i = 0; i < squarefree_factor_count; ++i)
    {
        if (integer_compare(output_stack, local_stack,
            integer_magnitude(local_stack, squarefree_factors[i]->coefficients[0]),
            integer_magnitude(local_stack, squarefree_factors[i]->coefficients
                [squarefree_factors[i]->coefficient_count - 1])) < 0)
        {
            array_reverse(squarefree_factors[i]->coefficients,
                squarefree_factors[i]->coefficient_count);
            struct IntegerPolynomial**reversed_irreducible_factors =
                stack_slot_allocate(local_stack,
                (squarefree_factor_count - 1) * sizeof(struct IntegerPolynomial*),
                    _Alignof(struct IntegerPolynomial*));
            size_t reversed_factor_count = squarefree_integer_polynomial_factor(local_stack,
                output_stack, squarefree_factors[i], reversed_irreducible_factors);
            for (size_t j = 0; j < reversed_factor_count; ++j)
            {
                array_reverse(reversed_irreducible_factors[j]->coefficients,
                    reversed_irreducible_factors[j]->coefficient_count);
                out[factor_count] =
                    integer_polynomial_copy(output_stack, reversed_irreducible_factors[j]);
                ++factor_count;
            }
        }
        else
        {
            factor_count += squarefree_integer_polynomial_factor(output_stack, local_stack,
                squarefree_factors[i], out);
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
}

struct IntegerPolynomial*integer_polynomial_derivative(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    return polynomial_derivative(integer_multiply, output_stack, local_stack,
        (struct Polynomial*)a);
}

struct RationalPolynomial*integer_polynomial_to_monic(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a)
{
    struct RationalPolynomial*out = stack_polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_reduced(output_stack, local_stack, a->coefficients[i],
            a->coefficients[a->coefficient_count - 1]);
    }
    return out;
}