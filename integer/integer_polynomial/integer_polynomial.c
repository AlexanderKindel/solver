#include "declarations.h"

struct IntegerPolynomial*integer_polynomial_copy(struct Stack*output_stack,
    struct IntegerPolynomial*a)
{
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_copy(output_stack, a->coefficients[i]);
    }
    return out;
}

void integer_polynomial_trim_leading_zeroes(struct IntegerPolynomial*a)
{
    for (size_t i = a->coefficient_count; i-- > 0;)
    {
        if (integer_equals(a->coefficients[i], &zero))
        {
            a->coefficient_count -= 1;
        }
        else
        {
            return;
        }
    }
}

bool integer_polynomial_equals(struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    if (a->coefficient_count != b->coefficient_count)
    {
        return false;
    }
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        if (!integer_equals(a->coefficients[i], b->coefficients[i]))
        {
            return false;
        }
    }
    return true;
}

struct IntegerPolynomial*integer_polynomial_add(struct Stack*restrict output_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    if (a->coefficient_count < b->coefficient_count)
    {
        SWAP(a, b, struct IntegerPolynomial*);
    }
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < b->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_add(output_stack, a->coefficients[i], b->coefficients[i]);
    }
    for (size_t i = b->coefficient_count; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_copy(output_stack, a->coefficients[i]);
    }
    integer_polynomial_trim_leading_zeroes(out);
    return out;
}

struct IntegerPolynomial*integer_polynomial_negate(struct Stack*restrict output_stack,
    struct IntegerPolynomial*a)
{
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_negate(output_stack, a->coefficients[i]);
    }
    return out;
}

struct IntegerPolynomial*integer_polynomial_subtract(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*minuend,
    struct IntegerPolynomial*subtrahend)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out = integer_polynomial_add(output_stack, minuend,
        integer_polynomial_negate(local_stack, subtrahend));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_integer_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Integer*b)
{
    if (integer_equals(b, &zero))
    {
        return polynomial_zero;
    }
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Integer*);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_multiply(output_stack, local_stack, a->coefficients[i], b);
    }
    return out;
}

struct RationalPolynomial*integer_polynomial_rational_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct Rational*b)
{
    if (b->numerator->sign == 0)
    {
        return polynomial_zero;
    }
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i] =
            rational_integer_multiply(output_stack, local_stack, b, a->coefficients[i]);
    }
    return out;
}

struct IntegerPolynomial*integer_polynomial_multiply(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial*b)
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
            out->coefficients[i + j] = integer_add(output_stack, out->coefficients[i + j],
                integer_multiply(local_stack, output_stack, a->coefficients[i],
                    b->coefficients[j]));
        }
    }
    integer_polynomial_trim_leading_zeroes(out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_integer_divide(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*dividend, struct Integer*divisor)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, dividend->coefficient_count, struct Integer*);
    for (size_t i = 0; i < dividend->coefficient_count; ++i)
    {
        out->coefficients[i] = integer_euclidean_divide(output_stack, local_stack,
            dividend->coefficients[i], divisor).quotient;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void integer_polynomial_copy_coefficients(struct Stack*output_stack, struct IntegerPolynomial*a)
{
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        a->coefficients[i] = integer_copy(output_stack, a->coefficients[i]);
    }
}

//Sets out->quotient and out->remainder to 0 when the coefficients of the quotient wouldn't all be
//integers.
struct IntegerPolynomialDivision integer_polynomial_euclidean_divide(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*dividend, struct IntegerPolynomial*divisor)
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
        dividend->coefficient_count * sizeof(struct Integer*));
    while (out.remainder->coefficient_count >= divisor->coefficient_count)
    {
        --out.remainder->coefficient_count;
        struct IntegerDivision division = integer_euclidean_divide(local_stack, output_stack,
            out.remainder->coefficients[out.remainder->coefficient_count],
            divisor->coefficients[divisor->coefficient_count - 1]);
        if (integer_equals(division.remainder, &zero))
        {
            out.quotient->coefficients[out.remainder->coefficient_count +
                1 - divisor->coefficient_count] = division.quotient;
            for (size_t i = 1; i < divisor->coefficient_count; ++i)
            {
                out.remainder->coefficients[out.remainder->coefficient_count - i] =
                    integer_subtract(local_stack, output_stack,
                        out.remainder->coefficients[out.remainder->coefficient_count - i],
                        integer_multiply(local_stack, output_stack, division.quotient,
                            divisor->coefficients[divisor->coefficient_count - i - 1]));
            }
        }
        else
        {
            out.quotient = 0;
            out.remainder = 0;
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
    }
    integer_polynomial_copy_coefficients(output_stack, out.quotient);
    integer_polynomial_trim_leading_zeroes(out.remainder);
    integer_polynomial_copy_coefficients(output_stack, out.remainder);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct IntegerPolynomial*integer_polynomial_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*base, struct Integer*exponent)
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
            out = integer_polynomial_multiply(local_stack, output_stack, out, base);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = integer_polynomial_copy(output_stack, out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = integer_polynomial_multiply(local_stack, output_stack, base, base);
    }
}

struct Integer*integer_polynomial_get_content(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return &one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*out = a->coefficients[0];
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        out = integer_get_gcd(local_stack, output_stack, out, a->coefficients[i]);
    }
    out = integer_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
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
    return (struct IntegerPolynomial*)integer_polynomial_get_extended_gcd(output_stack, local_stack,
        a, b).gcd;
}

struct PolynomialExtendedGCDInfo integer_polynomial_get_extended_gcd(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct IntegerPolynomial*a, struct IntegerPolynomial*b)
{
    if (b->coefficient_count > a->coefficient_count)
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct PolynomialExtendedGCDInfo out =
            integer_polynomial_get_extended_gcd(local_stack, output_stack, b, a);
        out.a_coefficient = 
            rational_polynomial_euclidean_divide(output_stack, local_stack,
                rational_polynomial_subtract(local_stack, output_stack,
                    integer_polynomial_to_rational_polynomial(local_stack,
                        (struct IntegerPolynomial*)out.gcd),
                    rational_polynomial_multiply(local_stack, output_stack,
                        integer_polynomial_to_rational_polynomial(local_stack, b),
                        out.a_coefficient)),
                integer_polynomial_to_rational_polynomial(local_stack, a)).quotient;
        out.gcd = (struct Polynomial*)integer_polynomial_copy(output_stack,
            (struct IntegerPolynomial*)out.gcd);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    if (b->coefficient_count == 0)
    {
        return (struct PolynomialExtendedGCDInfo) {
            (struct Polynomial*)integer_polynomial_copy(output_stack, a), rational_polynomial_one };
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*a_content = integer_polynomial_get_content(local_stack, output_stack, a);
    struct Integer*b_content = integer_polynomial_get_content(local_stack, output_stack, b);
    struct Integer*d = integer_get_gcd(local_stack, output_stack, a_content, b_content);
    struct Integer*g = &one;
    struct Integer*h = &one;
    struct IntegerPolynomial*previous_a_coefficient = integer_polynomial_one;
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
        struct IntegerPolynomialDivision division =
            integer_polynomial_euclidean_divide(local_stack, output_stack,
                integer_polynomial_integer_multiply(local_stack, output_stack,
                    previous_gcd_multiple, multiple),
                gcd_multiple);
        if (division.remainder->coefficient_count == 0)
        {
            struct Integer*gcd_content =
                integer_polynomial_get_content(local_stack, output_stack, gcd_multiple);
            struct Rational quotient = rational_reduce(local_stack, output_stack, d,
                integer_multiply(local_stack, output_stack, gcd_content, a_content));
            struct PolynomialExtendedGCDInfo out =
            { (struct Polynomial*)integer_polynomial_integer_multiply(output_stack, local_stack,
                integer_polynomial_integer_divide(local_stack, output_stack, gcd_multiple,
                    gcd_content), d),
                integer_polynomial_rational_multiply(output_stack, local_stack, a_coefficient,
                    &quotient) };
            local_stack->cursor = local_stack_savepoint;
            return out;
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
            h = integer_euclidean_divide(local_stack, output_stack,
                integer_exponentiate(local_stack, output_stack, g, degree),
                integer_exponentiate(local_stack, output_stack, h, h_exponent)).quotient;
        }
    }
}

struct IntegerPolynomial*integer_polynomial_get_derivative(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a)
{
    if (!a->coefficient_count)
    {
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct IntegerPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count - 1, struct Integer*);
    struct Integer*multiplier = &zero;
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        multiplier = integer_add(local_stack, multiplier, &one);
        out->coefficients[i - 1] =
            integer_multiply(output_stack, local_stack, a->coefficients[i], multiplier);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t integer_polynomial_squarefree_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial*a, struct IntegerPolynomial**out)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    struct IntegerPolynomial*b = a;
    struct IntegerPolynomial*c = integer_polynomial_get_derivative(local_stack, output_stack, a);
    a = integer_polynomial_get_gcd(local_stack, output_stack, b, c);
    struct IntegerPolynomialDivision division =
        integer_polynomial_euclidean_divide(local_stack, output_stack, b, a);
    b = division.quotient;
    do
    {
        division = integer_polynomial_euclidean_divide(local_stack, output_stack, c, a);
        c = integer_polynomial_subtract(local_stack, output_stack, division.quotient,
            integer_polynomial_get_derivative(local_stack, output_stack, b));
        a = integer_polynomial_get_gcd(output_stack, local_stack, b, c);
        if (a->coefficient_count > 1)
        {
            out[factor_count] = a;
            ++factor_count;
        }
        division = integer_polynomial_euclidean_divide(local_stack, output_stack, b, a);
        b = division.quotient;
    } while (b->coefficient_count > 1);
    local_stack->cursor = local_stack_savepoint;
    return factor_count;
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
            candidate->coefficients[i] = integer_euclidean_divide(local_stack, output_stack,
                candidate->coefficients[i], characteristic_power).remainder;
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
        struct IntegerPolynomialDivision division =
            integer_polynomial_euclidean_divide(local_stack, output_stack,
                integer_polynomial_integer_multiply(local_stack, output_stack, *a,
                    (*a)->coefficients[(*a)->coefficient_count - 1]), candidate);
        if (division.remainder && division.remainder->coefficient_count == 0)
        {
            factors[*factor_count] =
                integer_polynomial_get_primitive_part(output_stack, local_stack, candidate);
            while (true)
            {
                division = integer_polynomial_euclidean_divide(local_stack, output_stack, *a,
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
        modded_a_factors, modded_a, characteristic);
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
                integer_euclidean_divide(local_stack, output_stack, binomial_coefficient_numerator,
                    binomial_coefficient_denominator).quotient,
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
        product_of_unlifted_factors = modded_polynomial_euclidean_divide(local_stack, output_stack,
            product_of_unlifted_factors, modded_a_factors[i], characteristic).quotient;
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
            struct ExtendedGCDInfo gcd_info = modded_polynomial_get_extended_gcd(local_stack,
                output_stack, b_mod_prime, c_mod_prime, characteristic);
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
            struct IntegerPolynomialDivision division =
                modded_polynomial_euclidean_divide(local_stack, output_stack,
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
        SWAP(a->coefficients[a->coefficient_count - 1 - i], a->coefficients[i], struct Integer*);
    }
}

size_t primitive_integer_polynomial_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct IntegerPolynomial**out, struct IntegerPolynomial*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t factor_count = 0;
    if (!a->coefficients[0]->value_count)
    {
        out[0] = POLYNOMIAL_ALLOCATE(output_stack, 2, struct Integer*);
        out[0]->coefficients[0] = &zero;
        out[0]->coefficients[1] = &one;
        factor_count = 1;
        do
        {
            a = integer_polynomial_euclidean_divide(local_stack, output_stack, a, *out).quotient;
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
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
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
    struct RationalPolynomial*out =
        POLYNOMIAL_ALLOCATE(output_stack, a->coefficient_count, struct Rational);
    for (size_t i = 0; i < a->coefficient_count; ++i)
    {
        out->coefficients[i].numerator = integer_copy(output_stack, a->coefficients[i]);
        out->coefficients[i].denominator = &one;
    }
    return out;
}