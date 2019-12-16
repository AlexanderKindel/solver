﻿#include "declarations.h"

void*polynomial_allocate(struct Stack*output_stack, size_t coefficient_count)
{
    size_t*out = stack_slot_allocate(output_stack, sizeof(size_t) * (coefficient_count + 1),
        _Alignof(size_t));
    *out = coefficient_count;
    return out;
}

struct RationalPolynomial*integer_polynomial_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Rational*b)
{
    if (b->numerator->sign == 0)
    {
        return polynomial_zero;
    }
    struct RationalPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            rational_integer_multiply(output_stack, local_stack, b, a->coefficients[i]);
    }
    return out;
}

struct IntegerPolynomial*integer_polynomial_get_quotient(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor)
{
    struct IntegerPolynomialDivision division;
    integer_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return division.quotient;
}

struct IntegerPolynomial*integer_polynomial_get_remainder(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend,
    struct IntegerPolynomial*divisor)
{
    struct IntegerPolynomialDivision division;
    integer_polynomial_euclidean_divide(output_stack, local_stack, &division, dividend, divisor);
    return division.remainder;
}

struct IntegerPolynomial*integer_polynomial_get_primitive_part(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_integer_divide(output_stack, local_stack, a,
        integer_polynomial_get_content(local_stack, output_stack, a));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_get_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    struct PolynomialExtendedGCDInfo info;
    integer_polynomial_get_extended_gcd(output_stack, local_stack, &info, a, b);
    return (struct IntegerPolynomial*)info.gcd;
}

void integer_polynomial_get_extended_gcd(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct PolynomialExtendedGCDInfo*out,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    if (b->coefficient_count > a->coefficient_count)
    {
        void*local_stack_savepoint = local_stack->cursor;
        integer_polynomial_get_extended_gcd(local_stack, output_stack, out, b, a);
        out->a_coefficient = 
            rational_polynomial_get_quotient(output_stack, local_stack,
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
    struct Integer*a_content = integer_polynomial_get_content(local_stack, output_stack, a);
    struct Integer*b_content = integer_polynomial_get_content(local_stack, output_stack, b);
    struct Integer*d = integer_get_gcd(local_stack, output_stack, a_content, b_content);
    struct Integer*g = &one;
    struct Integer*h = &one;
    struct IntegerPolynomial*previous_a_coefficient = &integer_polynomial_one;
    struct IntegerPolynomial*previous_gcd_multiple =
        integer_polynomial_integer_divide(local_stack, output_stack, a, a_content);
    struct IntegerPolynomial*a_coefficient = polynomial_zero;
    struct IntegerPolynomial*gcd_multiple =
        integer_polynomial_integer_divide(local_stack, output_stack, b, b_content);
    while (true)
    {
        struct Integer*degree = size_t_to_integer(local_stack,
            previous_gcd_multiple->coefficient_count - gcd_multiple->coefficient_count);
        struct Integer*multiple = integer_exponentiate(local_stack, output_stack,
            gcd_multiple->coefficients[gcd_multiple->coefficient_count - 1],
            integer_add(local_stack, degree, &one));
        struct IntegerPolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, previous_gcd_multiple,
                multiple), gcd_multiple);
        if (division.remainder->coefficient_count == 0)
        {
            struct Integer*gcd_content =
                integer_polynomial_get_content(local_stack, output_stack, gcd_multiple);
            out->gcd =
                (struct Polynomial*)integer_polynomial_integer_multiply(output_stack, local_stack,
                    integer_polynomial_integer_divide(local_stack, output_stack, gcd_multiple,
                        gcd_content), d);
            out->a_coefficient = integer_polynomial_rational_multiply(output_stack, local_stack,
                a_coefficient, rational_reduce(local_stack, output_stack, d,
                    integer_multiply(local_stack, output_stack, gcd_content, a_content)));
            local_stack->cursor = local_stack_savepoint;
            return;
        }
        struct Integer*divisor = integer_multiply(local_stack, output_stack, g,
            integer_exponentiate(local_stack, output_stack, h, degree));
        struct IntegerPolynomial*next_a_coefficient =
            integer_polynomial_subtract(local_stack, output_stack,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    previous_a_coefficient, multiple),
                integer_polynomial_multiply(local_stack, output_stack, division.quotient,
                    a_coefficient));
        previous_a_coefficient = a_coefficient;
        previous_gcd_multiple = gcd_multiple;
        a_coefficient = integer_polynomial_integer_divide(local_stack, output_stack,
            next_a_coefficient, divisor);
        gcd_multiple = integer_polynomial_integer_divide(local_stack, output_stack,
            division.remainder, divisor);
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
            h = integer_get_quotient(local_stack, output_stack,
                integer_exponentiate(local_stack, output_stack, g, degree),
                integer_exponentiate(local_stack, output_stack, h, h_exponent));
        }
    }
}

bool factor_found(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
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
            candidate->coefficients[i] = integer_get_remainder(local_stack, output_stack,
                candidate->coefficients[i], characteristic_power);
            if (candidate->coefficients[i]->sign < 0)
            {
                candidate->coefficients[i] =
                    integer_add(local_stack, candidate->coefficients[i], characteristic_power);
            }
            if (integer_compare(output_stack, local_stack,
                integer_double(local_stack, candidate->coefficients[i]), characteristic_power) > 0)
            {
                candidate->coefficients[i] = integer_subtract(local_stack, output_stack,
                    candidate->coefficients[i], characteristic_power);
            }
        }
        struct IntegerPolynomialDivision division;
        integer_polynomial_euclidean_divide(local_stack, output_stack, &division,
            integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                (*a)->coefficients[(*a)->coefficient_count - 1]), candidate);
        if (division.remainder && division.remainder->coefficient_count == 0)
        {
            factors[*factor_count] =
                integer_polynomial_get_primitive_part(output_stack, local_stack, candidate);
            while (true)
            {
                integer_polynomial_euclidean_divide(local_stack, output_stack, &division, *a,
                    factors[*factor_count]);
                if (!division.remainder || division.remainder->coefficient_count != 0)
                {
                    *a = integer_polynomial_copy(output_stack, *a);
                    *factor_count += 1;
                    local_stack->cursor = local_stack_savepoint;
                    return true;
                }
                *a = division.quotient;
            }
        }
        local_stack->cursor = local_stack_savepoint;
    }
    else
    {
        for (size_t i = 0; i <= component_count - components_to_add; ++i)
        {
            void*local_stack_savepoint = local_stack->cursor;
            bool out = factor_found(output_stack, local_stack,
                integer_polynomial_multiply(local_stack, output_stack, candidate,
                    factor_components[i]),
                components_to_add - 1, factor_components + i + 1, component_count - i - 1, a,
                characteristic_power, factors, factor_count);
            local_stack->cursor = local_stack_savepoint;
            if (out)
            {
                factor_components[i] = factor_components[component_count - 1];
                return true;
            }
        }
    }
    return false;
}

size_t squarefree_integer_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*modded_a;
    struct IntegerPolynomial*gcd;
    next_prime = prime_stack.start;
    struct Integer*characteristic;
    while (true)
    {
        characteristic = get_next_prime(output_stack, local_stack);
        modded_a = modded_polynomial_reduce(local_stack, output_stack, a, characteristic);
        if (modded_a->coefficient_count == a->coefficient_count)
        {
            gcd = modded_polynomial_get_gcd(local_stack, output_stack, modded_a,
                modded_polynomial_reduce(local_stack, output_stack,
                    integer_polynomial_get_derivative(local_stack, output_stack, modded_a),
                    characteristic),
                characteristic);
            if (gcd->coefficient_count == 1)
            {
                break;
            }
        }
    }
    modded_a = modded_polynomial_get_monic(local_stack, output_stack, modded_a, characteristic);
    struct IntegerPolynomial**modded_a_factors = ARRAY_ALLOCATE(local_stack,
        modded_a->coefficient_count - 1, struct IntegerPolynomial*);
    size_t modded_a_factor_count = squarefree_modded_polynomial_factor(local_stack, output_stack,
        modded_a, characteristic, modded_a_factors);
    struct Integer*coefficient_bound = &zero;
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        coefficient_bound = integer_add(local_stack, coefficient_bound,
            integer_multiply(local_stack, output_stack, a->coefficients[i], a->coefficients[i]));
    }
    size_t factor_degree = (a->coefficient_count - 1) / 2;
    struct Integer*binomial_coefficient_numerator = &one;
    struct Integer*binomial_coefficient_denominator = &one;
    for (size_t i = (factor_degree - 1) / 2; i > 0; --i)
    {
        binomial_coefficient_numerator =
            integer_multiply(local_stack, output_stack, binomial_coefficient_numerator,
                size_t_to_integer(local_stack, factor_degree - i));
        binomial_coefficient_denominator = integer_multiply(local_stack, output_stack,
            binomial_coefficient_denominator, size_t_to_integer(local_stack, i));
    }
    struct Integer*leading_coefficient_magnitude =
        integer_get_magnitude(local_stack, a->coefficients[a->coefficient_count - 1]);
    coefficient_bound = integer_double(local_stack,
        integer_multiply(local_stack, output_stack, leading_coefficient_magnitude,
            integer_multiply(local_stack, output_stack,
                integer_get_quotient(local_stack, output_stack, binomial_coefficient_numerator,
                    binomial_coefficient_denominator),
                integer_add(local_stack, leading_coefficient_magnitude,
                    integer_take_square_root(local_stack, output_stack, coefficient_bound)))));
    struct Integer*e = &one;
    struct Integer*characteristic_power = characteristic;
    while (integer_compare(output_stack, local_stack, characteristic_power, coefficient_bound) < 0)
    {
        e = integer_add(local_stack, e, &one);
        characteristic_power =
            integer_multiply(local_stack, output_stack, characteristic_power, characteristic);
    }
    struct IntegerPolynomial*product_of_unlifted_factors = modded_a;
    struct IntegerPolynomial*unmodded_b_times_c = modded_polynomial_get_monic(local_stack,
        output_stack, modded_polynomial_reduce(local_stack, output_stack, a, characteristic_power),
        characteristic_power);
    struct IntegerPolynomial**lifted_factors =
        ARRAY_ALLOCATE(local_stack, modded_a_factor_count, struct IntegerPolynomial*);
    size_t lifted_factor_count = 0;
    modded_a_factor_count -= 1;
    for (size_t i = 0; i < modded_a_factor_count; ++i)
    {
        product_of_unlifted_factors = modded_polynomial_get_quotient(local_stack, output_stack,
            product_of_unlifted_factors, modded_a_factors[i], characteristic);
        struct Integer*lift_characteristic = characteristic;
        struct IntegerPolynomial*b = modded_a_factors[i];
        struct IntegerPolynomial*c = product_of_unlifted_factors;
        while (integer_compare(output_stack, local_stack, lift_characteristic,
            characteristic_power) < 0)
        {
            struct IntegerPolynomial*b_mod_prime =
                modded_polynomial_reduce(local_stack, output_stack, b, characteristic);
            struct IntegerPolynomial*c_mod_prime =
                modded_polynomial_reduce(local_stack, output_stack, c, characteristic);
            struct ExtendedGCDInfo gcd_info;
            modded_polynomial_get_extended_gcd(local_stack, output_stack, &gcd_info, b_mod_prime,
                c_mod_prime, characteristic);
            struct Integer*reciprocal = modded_integer_get_reciprocal(local_stack, output_stack,
                ((struct IntegerPolynomial*)gcd_info.gcd)->coefficients[0], characteristic);
            gcd_info.a_coefficient = modded_polynomial_modded_integer_multiply(local_stack,
                output_stack, (struct IntegerPolynomial*)gcd_info.a_coefficient, reciprocal,
                characteristic);
            gcd_info.b_coefficient = modded_polynomial_modded_integer_multiply(local_stack,
                output_stack, (struct IntegerPolynomial*)gcd_info.b_coefficient, reciprocal,
                characteristic);
            struct IntegerPolynomial*f = modded_polynomial_reduce(local_stack, output_stack,
                integer_polynomial_integer_divide(local_stack, output_stack,
                    integer_polynomial_subtract(local_stack, output_stack, unmodded_b_times_c,
                        integer_polynomial_multiply(local_stack, output_stack, b, c)),
                    lift_characteristic), characteristic);
            struct IntegerPolynomialDivision division;
            modded_polynomial_euclidean_divide(local_stack, output_stack, &division,
                modded_polynomial_multiply(local_stack, output_stack, f, gcd_info.b_coefficient,
                    characteristic),
                b_mod_prime, characteristic);
            b = integer_polynomial_add(local_stack, b,
                integer_polynomial_integer_multiply(local_stack, output_stack, division.remainder,
                    lift_characteristic));
            c = integer_polynomial_add(local_stack, c,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    modded_polynomial_add(local_stack, output_stack,
                        modded_polynomial_multiply(local_stack, output_stack,
                            (struct IntegerPolynomial*)gcd_info.a_coefficient, f, characteristic),
                        modded_polynomial_multiply(local_stack, output_stack, c_mod_prime,
                            division.quotient, characteristic),
                        characteristic),
                    lift_characteristic));
            lift_characteristic =
                integer_multiply(local_stack, output_stack, lift_characteristic, characteristic);
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
        for (size_t i = 0; i <= lifted_factor_count - combination_size;)
        {
            if (factor_found(output_stack, local_stack, lifted_factors[i], combination_size - 1,
                lifted_factors + i + 1, lifted_factor_count - i - 1, &a, characteristic_power, out,
                &factor_count))
            {
                lifted_factors[i] = lifted_factors[lifted_factor_count - 1];
                lifted_factor_count -= combination_size;
            }
            else
            {
                ++i;
            }
        }
        ++combination_size;
    }
    if (2 * combination_size == lifted_factor_count)
    {
        lifted_factor_count -= 1;
        factor_found(output_stack, local_stack, lifted_factors[lifted_factor_count],
            combination_size - 1, lifted_factors, lifted_factor_count, &a, characteristic_power,
            out, &factor_count);
    }
    if (a->coefficient_count > 1)
    {
        out[factor_count] = integer_polynomial_get_primitive_part(output_stack, local_stack, a);
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

void integer_polynomial_reverse(struct IntegerPolynomial*a)
{
    for (size_t i = 0; i < a->coefficient_count / 2; ++i)
    {
        POINTER_SWAP(a->coefficients[a->coefficient_count - 1 - i], a->coefficients[i]);
    }
}

size_t primitive_integer_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial**out)
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
            a = integer_polynomial_get_quotient(local_stack, output_stack, a, *out);
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
            integer_get_magnitude(local_stack, squarefree_factors[i]->coefficients[0]),
            integer_get_magnitude(local_stack, squarefree_factors[i]->coefficients
                [squarefree_factors[i]->coefficient_count - 1])) < 0)
        {
            integer_polynomial_reverse(squarefree_factors[i]);
            size_t reversed_factor_count = squarefree_integer_polynomial_factor(output_stack,
                local_stack, squarefree_factors[i], out + factor_count);
            for (size_t j = 0; j < reversed_factor_count; ++j)
            {
                integer_polynomial_reverse(out[factor_count]);
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

struct RationalPolynomial*integer_polynomial_get_monic(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a)
{
    struct RationalPolynomial*out = polynomial_allocate(output_stack, a->coefficient_count);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = rational_reduce(output_stack, local_stack, a->coefficients[i],
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