#include "declarations.h"

struct Number*number_surd_initialize(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*radicand, struct Integer*index)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '^';
    out->radicand = number_copy(output_stack, radicand);
    out->index = integer_copy(output_stack, index);
    size_t surd_index = integer_to_size_t(index);
    struct RationalPolynomial*annulling_polynomial = POLYNOMIAL_ALLOCATE(local_stack,
        surd_index * (radicand->minimal_polynomial->coefficient_count - 1) + 1, struct Rational);
    annulling_polynomial->coefficients[0] = radicand->minimal_polynomial->coefficients[0];
    for (size_t i = 1; i < radicand->minimal_polynomial->coefficient_count; ++i)
    {
        size_t annulling_polynomial_coefficient_index = i * surd_index;
        annulling_polynomial->coefficients[annulling_polynomial_coefficient_index] =
            radicand->minimal_polynomial->coefficients[i];
        for (size_t j = 1; j < surd_index; ++j)
        {
            annulling_polynomial->coefficients[annulling_polynomial_coefficient_index - j] =
                rational_zero;
        }
    }
    out->minimal_polynomial = number_annulling_polynomial_to_minimal_polynomial(output_stack,
        local_stack, annulling_polynomial, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Rational number_get_rational_factor(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        return rational_get_magnitude(output_stack, &a->value);
    case '^':
        return rational_one;
    case '*':
        if (a->elements[0]->operation == 'r')
        {
            return rational_get_magnitude(output_stack, &a->elements[0]->value);
        }
        return rational_one;
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational out = number_get_rational_factor(local_stack, output_stack, a->elements[0]);
        for (size_t i = 1; i < a->element_count; ++i)
        {
            struct Rational term_factor =
                number_get_rational_factor(local_stack, output_stack, a->elements[i]);
            out = rational_reduce(local_stack, output_stack,
                integer_get_gcd(local_stack, output_stack, out.numerator, term_factor.numerator),
                integer_get_lcm(local_stack, output_stack, out.denominator,
                    term_factor.denominator));
        }
        out = rational_copy(output_stack, &out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    }
    crash("Number operation not recognized.");
}

struct Number*number_integer_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*exponent)
{
    if (!exponent->value_count)
    {
        return &number_one;
    }
    if (exponent->sign < 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        base = number_get_reciprocal(local_stack, output_stack, base);
        if (!base)
        {
            local_stack->cursor = local_stack_savepoint;
            return 0;
        }
        struct Number*out = number_integer_exponentiate(output_stack, local_stack, base,
            integer_negate(local_stack, exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    switch (base->operation)
    {
    case 'r':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational value_exponentiation =
            rational_exponentiate(local_stack, output_stack, base->value, exponent);
        struct Number*out = number_rational_initialize(output_stack, &value_exponentiation);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '^':
        return number_rational_exponentiate(output_stack, local_stack, base->radicand,
            &(struct Rational){ exponent, base->index });
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = &number_one;
        for (size_t i = base->element_count; i-- > 1;)
        {
            out = number_multiply(local_stack, output_stack, out,
                number_integer_exponentiate(local_stack, output_stack, base->elements[i],
                    exponent));
        }
        out = number_multiply(output_stack, local_stack, out,
            number_integer_exponentiate(local_stack, output_stack, base->elements[0], exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
        return number_exponentiate(output_stack, local_stack, base, exponent);
    }
    crash("Number operation not recognized.");
}

struct Number*number_rational_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Rational*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out =
        number_integer_exponentiate(local_stack, output_stack, base, exponent->numerator);
    if (out)
    {
        out = number_take_root(output_stack, local_stack, out, exponent->denominator);
    }
    else
    {
        out = 0;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_exponentiate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*exponent)
{
    if (!exponent->value_count)
    {
        return &number_one;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = &number_one;
    while (true)
    {
        if (exponent->value[0] & 1)
        {
            out = number_multiply(local_stack, output_stack, out, base);
        }
        exponent = integer_halve(local_stack, exponent);
        if (!exponent->value_count)
        {
            out = number_copy(output_stack, out);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        base = number_multiply(local_stack, output_stack, base, base);
    }
}

struct Number*number_take_root(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*base, struct Integer*index)
{
    if (!index->value_count)
    {
        return &number_one;
    }
    if (integer_equals(index, &one))
    {
        return number_copy(output_stack, base);
    }
    void*local_stack_savepoint = local_stack->cursor;
    base = number_eliminate_linear_dependencies(local_stack, output_stack, base);
    switch (base->operation)
    {
    case 'r':
    {
        if (!integer_equals(base->value.denominator, &one))
        {
            struct Number*out = number_rational_multiply(output_stack, local_stack,
                number_take_root(local_stack, output_stack,
                    number_rational_initialize(local_stack,
                        &(struct Rational) { integer_multiply(local_stack, output_stack,
                            base->value.numerator,
                            integer_exponentiate(local_stack, output_stack,
                                base->value.denominator,
                                integer_add(local_stack, index, INT(1, -1)))), &one}), index),
                &(struct Rational) { &one, base->value.denominator });
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Integer*radicand = integer_get_magnitude(local_stack, base->value.numerator);
        struct Factor*factors;
        size_t factor_count = integer_factor(local_stack, output_stack, &factors, radicand);
        struct Integer*coefficient = &one;
        struct Integer*multiplicity_gcd = index;
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct IntegerDivision multiplicity_reduction = integer_euclidean_divide(local_stack,
                output_stack, factors[factor_index].multiplicity, index);
            factors[factor_index].multiplicity = multiplicity_reduction.remainder;
            coefficient = integer_multiply(local_stack, output_stack, coefficient,
                integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                    multiplicity_reduction.quotient));
            multiplicity_gcd = integer_get_gcd(local_stack, output_stack, multiplicity_gcd,
                factors[factor_index].multiplicity);
        }
        struct Integer*new_index;
        if (base->value.numerator->sign > 0)
        {
            new_index = integer_euclidean_divide(local_stack, output_stack, index,
                multiplicity_gcd).quotient;
            if (integer_equals(new_index, &one))
            {
                struct Rational*rational_radicand = ALLOCATE(output_stack, struct Rational);
                rational_radicand->numerator = integer_copy(output_stack, coefficient);
                rational_radicand->denominator = &one;
                struct Number*out = number_rational_initialize(output_stack, rational_radicand);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
        }
        else
        {
            new_index = index;
        }
        radicand = integer_initialize(local_stack, 1, base->value.numerator->sign);
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_divide(local_stack,
                output_stack, factors[factor_index].multiplicity, multiplicity_gcd).quotient;
            struct Integer*exponentiation = integer_exponentiate(local_stack, output_stack,
                factors[factor_index].value, reduced_multiplicity);
            radicand = integer_multiply(local_stack, output_stack, radicand, exponentiation);
        }
        struct Number*out = number_rational_multiply(output_stack, local_stack,
            number_surd_initialize(local_stack, output_stack,
                number_rational_initialize(local_stack, &(struct Rational){radicand, &one}),
                new_index),
            &(struct Rational){coefficient, &one});
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '^':
    {
        struct Number*out = number_take_root(output_stack, local_stack, base->radicand,
            integer_multiply(local_stack, output_stack, index, base->index));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        struct Number*out = &number_one;
        for (size_t i = base->element_count; i-- > 1;)
        {
            out = number_multiply(local_stack, output_stack, out,
                number_take_root(local_stack, output_stack, base->elements[i], index));
        }
        out = number_multiply(output_stack, local_stack, out,
            number_take_root(local_stack, output_stack, base->elements[0], index));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        struct Rational base_rational_factor =
            number_get_rational_factor(local_stack, output_stack, base);
        struct Factor*factors;
        size_t factor_count =
            integer_factor(local_stack, output_stack, &factors, base_rational_factor.numerator);
        struct Rational base_cancelling_rational_factor = { integer_exponentiate(local_stack,
            output_stack, base_rational_factor.denominator, index), &one };
        struct Rational product_rational_factor = { &one, base_rational_factor.denominator };
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_divide(local_stack,
                output_stack, factors[factor_index].multiplicity, index).quotient;
            base_cancelling_rational_factor.denominator =
                integer_multiply(local_stack, output_stack,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        integer_multiply(local_stack, output_stack, reduced_multiplicity, index)),
                    base_cancelling_rational_factor.denominator);
            product_rational_factor.numerator =
                integer_multiply(local_stack, output_stack, product_rational_factor.numerator,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        reduced_multiplicity));
        }
        struct Number*out = number_rational_multiply(output_stack, local_stack,
            number_surd_initialize(local_stack, output_stack,
                number_rational_multiply(local_stack, output_stack, base,
                    &base_cancelling_rational_factor), index), &product_rational_factor);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    }
    crash("Number operation not recognized.");
}