#include "declarations.h"

struct IntegerPolynomial*integer_polynomial_copy(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    return polynomial_copy(integer_copy, output_stack, (struct Polynomial*)a);
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

struct RationalPolynomial*integer_polynomial_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct IntegerPolynomial*a, struct Rational*b)
{
    if (b->numerator->sign == 0)
    {
        return (struct RationalPolynomial*)&polynomial_zero;
    }
    struct RationalPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            rational_integer_multiply(output_stack, local_stack, b, a->coefficients[i]);
    }
    return out;
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
        void*local_stack_savepoint = local_stack->cursor;
        integer_polynomial_extended_gcd(local_stack, output_stack, out, b, a);
        out->a_coefficient = 
            rational_polynomial_euclidean_quotient(output_stack, local_stack,
                rational_polynomial_subtract(local_stack, output_stack,
                    integer_polynomial_to_rational_polynomial(local_stack,
                        (struct IntegerPolynomial*)out->gcd),
                    rational_polynomial_multiply(local_stack, output_stack,
                        integer_polynomial_to_rational_polynomial(local_stack, b),
                        out->a_coefficient)),
                integer_polynomial_to_rational_polynomial(local_stack, a));
        out->gcd = (struct Polynomial*)integer_polynomial_copy(output_stack,
            (struct IntegerPolynomial*)out->gcd);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    if (b->coefficient_count == 0)
    {
        out->gcd = (struct Polynomial*)integer_polynomial_copy(output_stack, a);
        out->a_coefficient = &rational_polynomial_one;
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
    struct IntegerPolynomial*a_coefficient = (struct IntegerPolynomial*)&polynomial_zero;
    struct IntegerPolynomial*gcd_multiple =
        integer_polynomial_integer_divide(local_stack, output_stack, b, b_content);
    while (true)
    {
        struct Integer*degree = integer_from_size_t(local_stack,
            previous_gcd_multiple->coefficient_count - gcd_multiple->coefficient_count);
        struct Integer*multiple = integer_exponentiate(local_stack, output_stack,
            gcd_multiple->coefficients[gcd_multiple->coefficient_count - 1],
            integer_add(local_stack, degree, &one));
        struct PolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, previous_gcd_multiple,
                multiple), gcd_multiple);
        if (division.remainder->coefficient_count == 0)
        {
            struct Integer*gcd_content =
                integer_polynomial_content(local_stack, output_stack, gcd_multiple);
            out->gcd =
                (struct Polynomial*)integer_polynomial_integer_multiply(output_stack, local_stack,
                    integer_polynomial_integer_divide(local_stack, output_stack, gcd_multiple,
                        gcd_content), d);
            struct Rational*multiple = rational_reduced(local_stack, output_stack, d, gcd_content);
            out->a_coefficient = integer_polynomial_rational_multiply(output_stack, local_stack,
                a_coefficient, multiple);
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        struct Integer*divisor = integer_multiply(local_stack, output_stack, g,
            integer_exponentiate(local_stack, output_stack, h, degree));
        struct IntegerPolynomial*next_a_coefficient =
            integer_polynomial_subtract(local_stack, output_stack,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    previous_a_coefficient, multiple),
                integer_polynomial_multiply(local_stack, output_stack,
                    (struct IntegerPolynomial*)division.quotient, a_coefficient));
        previous_a_coefficient = a_coefficient;
        previous_gcd_multiple = gcd_multiple;
        a_coefficient = integer_polynomial_integer_divide(local_stack, output_stack,
            next_a_coefficient, divisor);
        gcd_multiple = integer_polynomial_integer_divide(local_stack, output_stack,
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

bool find_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct IntegerPolynomial*candidate, size_t components_to_add,
    struct IntegerPolynomial**factor_components, size_t component_count,
    struct IntegerPolynomial**a, struct Integer*characteristic_power,
    struct IntegerPolynomial**factors, size_t*factor_count)
{
    if (components_to_add == 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        candidate = integer_polynomial_integer_multiply(local_stack, output_stack, candidate,
            (*a)->coefficients[(*a)->coefficient_count - 1]);
        for (size_t i = 0; i < candidate->coefficient_count; ++i)
        {
            candidate->coefficients[i] = integer_euclidean_remainder(local_stack, output_stack,
                candidate->coefficients[i], characteristic_power);
            if (candidate->coefficients[i]->sign < 0)
            {
                candidate->coefficients[i] =
                    integer_add(local_stack, candidate->coefficients[i], characteristic_power);
            }
            if (integer_compare(output_stack, local_stack,
                integer_doubled(local_stack, candidate->coefficients[i]), characteristic_power) > 0)
            {
                candidate->coefficients[i] = integer_subtract(local_stack, output_stack,
                    candidate->coefficients[i], characteristic_power);
            }
        }
        struct PolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, *a,
            (*a)->coefficients[(*a)->coefficient_count - 1]), candidate);
        if (division.remainder && division.remainder->coefficient_count == 0)
        {
            factors[*factor_count] =
                integer_polynomial_primitive_part(output_stack, local_stack, candidate);
            while (true)
            {
                integer_polynomial_euclidean_divide(local_stack, output_stack, &division, *a,
                    factors[*factor_count]);
                if (!division.remainder || division.remainder->coefficient_count != 0)
                {
                    *factor_count += 1;
                    local_stack->cursor = local_stack_savepoint;
                    return true;
                }
                *a = (struct IntegerPolynomial*)division.quotient;
            }
        }
        local_stack->cursor = local_stack_savepoint;
    }
    else
    {
        for (size_t i = 0; i <= component_count - components_to_add; ++i)
        {
            void*local_stack_savepoint = local_stack->cursor;
            bool factor_found = find_factor(output_stack, local_stack,
                integer_polynomial_multiply(local_stack, output_stack, candidate,
                    factor_components[i]), components_to_add - 1, factor_components + 1,
                component_count - 1, a, characteristic_power, factors, factor_count);
            local_stack->cursor = local_stack_savepoint;
            if (factor_found)
            {
                factor_components[i] = factor_components[component_count - 1];
                return factor_found;
            }
        }
    }
    return false;
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
    struct IntegerPolynomial**modded_a_factors = ARRAY_ALLOCATE(local_stack,
        modded_a->coefficient_count - 1, struct IntegerPolynomial*);
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
    struct IntegerPolynomial**lifted_factors =
        ARRAY_ALLOCATE(local_stack, modded_a_factor_count, struct IntegerPolynomial*);
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
    size_t factor_count = 0;
    size_t combination_size = 1;
    while (2 * combination_size < lifted_factor_count)
    {
        while (true)
        {
            bool factor_found = find_factor(output_stack, local_stack, &integer_polynomial_one,
                combination_size, lifted_factors, lifted_factor_count, &a, prime_power, out,
                &factor_count);
            if (!factor_found)
            {
                break;
            }
            lifted_factor_count -= combination_size;
        }
        ++combination_size;
    }
    if (2 * combination_size == lifted_factor_count)
    {
        find_factor(output_stack, local_stack, lifted_factors[lifted_factor_count - 1],
            combination_size - 1, lifted_factors, lifted_factor_count - 1, &a, prime_power, out,
            &factor_count);
    }
    if (a->coefficient_count > 1)
    {
        out[factor_count] = integer_polynomial_primitive_part(output_stack, local_stack, a);
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
        out[0] = polynomial_allocate(output_stack, 2);
        out[0]->coefficients[0] = &zero;
        out[0]->coefficients[1] = &one;
        factor_count = 1;
        do
        {
            a = integer_polynomial_euclidean_quotient(local_stack, output_stack, a, *out);
        } while (!a->coefficients[0]->value_count);
    }
    if (a->coefficient_count < 2)
    {
        local_stack->cursor = local_stack_savepoint;
        return factor_count;
    }
    struct IntegerPolynomial**squarefree_factors = ARRAY_ALLOCATE(local_stack,
        a->coefficient_count - 1, struct IntegerPolynomial*);
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
            size_t reversed_factor_count = squarefree_integer_polynomial_factor(output_stack,
                local_stack, squarefree_factors[i], out + factor_count);
            for (size_t j = 0; j < reversed_factor_count; ++j)
            {
                array_reverse(out[factor_count]->coefficients,
                    out[factor_count]->coefficient_count);
                ++factor_count;
            }
        }
        else
        {
            factor_count += squarefree_integer_polynomial_factor(output_stack, local_stack,
                squarefree_factors[i], out + factor_count);
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
    struct RationalPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_reduced(output_stack, local_stack, a->coefficients[i],
            a->coefficients[a->coefficient_count - 1]);
    }
    return out;
}

struct RationalPolynomial*integer_polynomial_to_rational_polynomial(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    struct RationalPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = ALLOCATE(output_stack, struct Rational);
        out->coefficients[i]->numerator = integer_copy(output_stack, a->coefficients[i]);
        out->coefficients[i]->denominator = &one;
    }
    return out;
}