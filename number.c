#include "declarations.h"

struct Number*number_copy(struct Stack*output_stack, struct Number*a)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = a->operation;
    if (a->minimal_polynomial)
    {
        out->minimal_polynomial = rational_polynomial_copy(output_stack, a->minimal_polynomial);
    }
    else
    {
        out->minimal_polynomial = 0;
    }
    switch (a->operation)
    {
    case 'r':
        out->value = rational_copy(output_stack, a->value);
        break;
    case '^':
        out->radicand = number_copy(output_stack, a->radicand);
        out->index = integer_copy(output_stack, a->index);
        break;
    case '+':
        if (a->generator)
        {
            out->generator = number_copy(output_stack, a->generator);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, a->element_count, struct RationalPolynomial*);
            for (size_t i = 0; i < a->element_count; ++i)
            {
                out->terms_in_terms_of_generator[i] =
                    rational_polynomial_copy(output_stack, a->terms_in_terms_of_generator[i]);
            }
        }
        else
        {
            out->generator = 0;
        }
    case '*':
        out->element_count = a->element_count;
        out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
        for (size_t i = 0; i < a->element_count; ++i)
        {
            out->elements[i] = number_copy(output_stack, a->elements[i]);
        }
    }
    return out;
}

struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = 'r';
    out->value = rational_copy(output_stack, value);
    out->minimal_polynomial = polynomial_allocate(output_stack, 2);
    out->minimal_polynomial->coefficients[0] = ALLOCATE(output_stack, struct Rational);
    out->minimal_polynomial->coefficients[0]->numerator =
        integer_negative(output_stack, value->numerator);
    out->minimal_polynomial->coefficients[0]->denominator =
        integer_copy(output_stack, value->denominator);
    out->minimal_polynomial->coefficients[1] = &rational_one;
    return out;
}

struct Number*number_surd_initialize(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*radicand, struct Integer*index)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '^';
    out->radicand = number_copy(output_stack, radicand);
    out->index = integer_copy(output_stack, index);
    size_t surd_index = integer_to_size_t(index);
    struct RationalPolynomial*annulling_polynomial = polynomial_allocate(local_stack,
        surd_index * (radicand->minimal_polynomial->coefficient_count - 1) + 1);
    annulling_polynomial->coefficients[0] = radicand->minimal_polynomial->coefficients[0];
    for (size_t i = 1; i < radicand->minimal_polynomial->coefficient_count; ++i)
    {
        size_t annulling_polynomial_coefficient_index = i * surd_index;
        annulling_polynomial->coefficients[annulling_polynomial_coefficient_index] =
            radicand->minimal_polynomial->coefficients[i];
        for (size_t j = 1; j < surd_index; ++j)
        {
            annulling_polynomial->coefficients[annulling_polynomial_coefficient_index - j] =
                &rational_zero;
        }
    }
    out->minimal_polynomial = number_minimal_polynomial_from_annulling_polynomial(output_stack,
        local_stack, annulling_polynomial, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

int8_t sum_or_product_formal_compare(struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct Number*b)
{
    if (a->element_count > b->element_count)
    {
        return 1;
    }
    if (a->element_count < b->element_count)
    {
        return -1;
    }
    for (size_t i = 0; i < a->element_count; ++i)
    {
        int8_t out = number_formal_compare(stack_a, stack_b, a->elements[i], b->elements[i]);
        if (out)
        {
            return out;
        }
    }
    return 0;
}

int8_t number_formal_compare(struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (b->operation == 'r')
        {
            return rational_compare(stack_a, stack_b, a->value, b->value);
        }
        return -1;
    case '^':
        switch (b->operation)
        {
        case 'r':
            return 1;
        case '^':
        {
            int8_t out = integer_compare(stack_a, stack_b, a->index, b->index);
            if (out)
            {
                return -out;
            }
            return number_formal_compare(stack_a, stack_b, a->radicand, b->radicand);
        }
        }
        return -1;
    case '*':
        switch (b->operation)
        {
        case '+':
            return -1;
        case '*':
            return sum_or_product_formal_compare(stack_a, stack_b, a, b);
        }
        return 1;
    case '+':
        if (b->operation != '+')
        {
            return 1;
        }
        return sum_or_product_formal_compare(stack_a, stack_b, a, b);
    }
    crash("Number operation not recognized.");
}

struct RationalPolynomial*number_minimal_polynomial_from_annulling_polynomial(
    struct Stack*output_stack, struct Stack*local_stack,
    struct RationalPolynomial*annulling_polynomial, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial**candidates = ARRAY_ALLOCATE(local_stack,
        annulling_polynomial->coefficient_count, struct RationalPolynomial*);
    size_t candidate_count =
        rational_polynomial_factor(local_stack, output_stack, annulling_polynomial, candidates);
    if (candidate_count == 1)
    {
        struct RationalPolynomial*out = rational_polynomial_copy(output_stack, candidates[0]);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct FloatInterval zero_estimate = { &float_zero, &float_zero };
    struct RectangularEstimate a_estimate;
    number_rectangular_estimate(local_stack, output_stack, &a_estimate, a, &rational_one);
    struct Rational interval_size_for_evaluation =
    { &one, integer_from_size_t(local_stack, annulling_polynomial->coefficient_count) };
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct RectangularEstimate evaluation_estimate;
            rational_polynomial_evaluate_at_rectangular_estimate(local_stack, output_stack,
                &evaluation_estimate, candidates[i], &a_estimate, &interval_size_for_evaluation);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &evaluation_estimate,
                &(struct RectangularEstimate){ &zero_estimate, &zero_estimate }))
            {
                candidate_count -= 1;
                candidates[i] = candidates[candidate_count];
                if (candidate_count == 1)
                {
                    struct RationalPolynomial*out =
                        rational_polynomial_copy(output_stack, candidates[0]);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
            else
            {
                ++i;
            }
        }
        number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &a_estimate, a);
        interval_size_for_evaluation.denominator =
            integer_doubled(local_stack, interval_size_for_evaluation.denominator);
    }
}

struct Number*number_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out;
    switch (a->operation)
    {
    case 'r':
        if (!a->value->numerator->value_count)
        {
            puts("Tried to divide by 0.");
            return 0;
        }
        out = number_rational_initialize(output_stack, rational_reciprocal(local_stack, a->value));
        break;
    case '^':
        out = number_divide(output_stack, local_stack,
            number_root(local_stack, output_stack,
                number_integer_exponentiate(local_stack, output_stack, a->radicand,
                    integer_add(local_stack, a->index, &INT(1, -))), a->index), a->radicand);
        break;
    case '*':
    {
        out = &number_one;
        for (size_t i = a->element_count; i-- > 1;)
        {
            out = number_multiply(local_stack, output_stack, out,
                number_reciprocal(local_stack, output_stack, a->elements[i]));
        }
        out = number_multiply(output_stack, local_stack, out,
            number_reciprocal(local_stack, output_stack, a->elements[0]));
        break;
    }
    case '+':
    {
        a = number_eliminate_linear_dependencies(local_stack, output_stack, a);
        struct Number**conjugates = number_conjugates(local_stack, output_stack, a);
        out = &number_one;
        for (size_t i = 1; i < a->minimal_polynomial->coefficient_count - 1; ++i)
        {
            out = number_multiply(local_stack, output_stack, out, conjugates[i]);
        }
        struct Rational*coefficient =
            rational_reciprocal(local_stack, a->minimal_polynomial->coefficients[0]);
        if (a->minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient->numerator->sign *= -1;
        }
        out = number_rational_multiply(output_stack, local_stack, out, coefficient);
        break;
    }
    default:
        crash("Number operation not recognized.");
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_divide(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*dividend, struct Number*divisor)
{
    struct Number*reciprocal = number_reciprocal(local_stack, output_stack, divisor);
    if (!reciprocal)
    {
        return 0;
    }
    return number_multiply(output_stack, local_stack, dividend, reciprocal);
}

struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        return rational_magnitude(output_stack, a->value);
    case '^':
        return &rational_one;
    case '*':
        if (a->elements[0]->operation == 'r')
        {
            return rational_magnitude(output_stack, a->elements[0]->value);
        }
        return &rational_one;
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational*out = number_rational_factor(local_stack, output_stack, a->elements[0]);
        for (size_t i = 1; i < a->element_count; ++i)
        {
            struct Rational*term_factor =
                number_rational_factor(local_stack, output_stack, a->elements[i]);
            out = rational_reduced(local_stack, output_stack,
                integer_gcd(local_stack, output_stack, out->numerator, term_factor->numerator),
                integer_lcm(local_stack, output_stack, out->denominator, term_factor->denominator));
        }
        out = rational_copy(output_stack, out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_integer_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*base, struct Integer*exponent)
{
    if (!exponent->value_count)
    {
        return &number_one;
    }
    if (exponent->sign < 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        base = number_reciprocal(local_stack, output_stack, base);
        if (!base)
        {
            local_stack->cursor = local_stack_savepoint;
            return 0;
        }
        struct Number*out = number_integer_exponentiate(output_stack, local_stack, base,
            integer_negative(local_stack, exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    switch (base->operation)
    {
    case 'r':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_rational_initialize(output_stack,
            rational_exponentiate(local_stack, output_stack, base->value, exponent));
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
        return generic_exponentiate(&(struct RingOperations){number_copy, 0, 0, &number_one, 0,
            0, number_generic_multiply}, output_stack, local_stack, base, exponent, 0);
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_rational_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*base, struct Rational*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out =
        number_integer_exponentiate(local_stack, output_stack, base, exponent->numerator);
    if (out)
    {
        out = number_root(output_stack, local_stack, out, exponent->denominator);
    }
    else
    {
        out = 0;
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_root(struct Stack*output_stack, struct Stack*local_stack, struct Number*base,
    struct Integer*index)
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
        if (!integer_equals(base->value->denominator, &one))
        {
            struct Number*out = number_rational_multiply(output_stack, local_stack,
                number_root(local_stack, output_stack,
                    number_rational_initialize(local_stack,
                        &(struct Rational){integer_multiply(local_stack, output_stack,
                            base->value->numerator,
                            integer_exponentiate(local_stack, output_stack,
                                base->value->denominator,
                                integer_add(local_stack, index, &INT(1, -)))), &one}), index),
                &(struct Rational){&one, base->value->denominator});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Integer*radicand = integer_magnitude(local_stack, base->value->numerator);
        struct Factor*factors;
        size_t factor_count = integer_factor(local_stack, output_stack, &factors, radicand);
        struct Integer*coefficient = &one;
        struct Integer*multiplicity_gcd = index;
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct IntegerDivision multiplicity_reduction;
            integer_euclidean_divide(local_stack, output_stack, &multiplicity_reduction,
                factors[factor_index].multiplicity, index);
            factors[factor_index].multiplicity = multiplicity_reduction.remainder;
            coefficient = integer_multiply(local_stack, output_stack, coefficient,
                integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                    multiplicity_reduction.quotient));
            multiplicity_gcd = integer_gcd(local_stack, output_stack, multiplicity_gcd,
                factors[factor_index].multiplicity);
        }
        struct Integer*new_index;
        if (base->value->numerator->sign > 0)
        {
            new_index =
                integer_euclidean_quotient(local_stack, output_stack, index, multiplicity_gcd);
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
        radicand = integer_initialize(local_stack, 1, base->value->numerator->sign);
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_quotient(local_stack,
                output_stack, factors[factor_index].multiplicity, multiplicity_gcd);
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
        struct Number*out = number_root(output_stack, local_stack, base->radicand,
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
                number_root(local_stack, output_stack, base->elements[i], index));
        }
        out = number_multiply(output_stack, local_stack, out,
            number_root(local_stack, output_stack, base->elements[0], index));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        struct Rational*base_rational_factor =
            number_rational_factor(local_stack, output_stack, base);
        struct Factor*factors;
        size_t factor_count =
            integer_factor(local_stack, output_stack, &factors, base_rational_factor->numerator);
        struct Rational base_cancelling_rational_factor = { integer_exponentiate(local_stack,
            output_stack, base_rational_factor->denominator, index), &one };
        struct Rational product_rational_factor = { &one, base_rational_factor->denominator };
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_quotient(local_stack,
                output_stack, factors[factor_index].multiplicity, index);
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
    default:
        crash("Number operation not recognized.");
    }
}

struct Number**binary_sum_or_product_conjugates(struct Number*(operation)(struct Stack*,
    struct Stack*, struct Number*, struct Number*), struct Stack*output_stack,
    struct Stack*local_stack, struct Number**a_conjugates, size_t a_conjugate_count,
    struct Number**b_conjugates, size_t b_conjugate_count,
    struct RationalPolynomial*minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number**out =
        ARRAY_ALLOCATE(output_stack, minimal_polynomial->coefficient_count - 1, struct Number*);
    size_t conjugate_count = 0;
    for (size_t i = 0; i < a_conjugate_count; ++i)
    {
        for (size_t j = 0; j < b_conjugate_count; ++j)
        {
            void*output_stack_savepoint = output_stack->cursor;
            out[conjugate_count] = number_eliminate_linear_dependencies(output_stack, local_stack,
                operation(local_stack, output_stack, a_conjugates[i], b_conjugates[j]));
            if (rational_polynomial_equals(minimal_polynomial,
                out[conjugate_count]->minimal_polynomial))
            {
                ++conjugate_count;
                if (conjugate_count == minimal_polynomial->coefficient_count - 1)
                {
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
            else
            {
                output_stack->cursor = output_stack_savepoint;
            }
        }
    }
    crash("Not enough conjugates found.");
}

struct Number**sum_or_product_conjugates(struct Number*(operation)(struct Stack*, struct Stack*,
    struct Number*, struct Number*),
    struct RationalPolynomial*(get_minimal_polynomial)(struct Stack*, struct Stack*, struct Number*,
        struct RationalPolynomial*, struct RationalPolynomial*),
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number**out = number_conjugates(local_stack, output_stack, a->elements[0]);
    struct RationalPolynomial*left_minimal_polynomial = a->elements[0]->minimal_polynomial;
    size_t original_element_count = a->element_count;
    a->element_count = 2;
    while (a->element_count < original_element_count)
    {
        struct Number*right_element = a->elements[a->element_count - 1];
        struct RationalPolynomial*new_left_minimal_polynomial =
            get_minimal_polynomial(local_stack, output_stack, a, left_minimal_polynomial,
                right_element->minimal_polynomial);
        out = binary_sum_or_product_conjugates(operation, local_stack, output_stack, out,
            left_minimal_polynomial->coefficient_count - 1,
            number_conjugates(local_stack, output_stack, right_element),
            right_element->minimal_polynomial->coefficient_count - 1, new_left_minimal_polynomial);
        left_minimal_polynomial = new_left_minimal_polynomial;
        ++a->element_count;
    }
    out = binary_sum_or_product_conjugates(operation, output_stack, local_stack, out,
        left_minimal_polynomial->coefficient_count - 1,
        number_conjugates(local_stack, output_stack, a->elements[a->element_count - 1]),
        a->elements[a->element_count - 1]->minimal_polynomial->coefficient_count - 1,
        a->minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number**number_conjugates(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        struct Number**out = ALLOCATE(output_stack, struct Number*);
        out[0] = number_copy(output_stack, a);
        return out;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number**radicand_conjugates =
            number_conjugates(local_stack, output_stack, a->radicand);
        for (size_t i = 0; i < a->radicand->minimal_polynomial->coefficient_count - 1; ++i)
        {
            radicand_conjugates[i] =
                number_root(local_stack, output_stack, radicand_conjugates[i], a->index);
        }
        size_t roots_of_unity_count = integer_to_size_t(a->index);
        struct Number**roots_of_unity = get_roots_of_unity(output_stack, local_stack, a->index);
        struct Number**out =
            binary_sum_or_product_conjugates(number_multiply, output_stack, local_stack,
                radicand_conjugates, a->radicand->minimal_polynomial->coefficient_count - 1,
                roots_of_unity, roots_of_unity_count, a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
        return sum_or_product_conjugates(number_multiply, product_minimal_polynomial, output_stack,
            local_stack, a);
    case '+':
        return sum_or_product_conjugates(number_add, sum_minimal_polynomial, output_stack,
            local_stack, a);
    default:
        crash("Number operation not recognized.");
    }
}

size_t trailing_factor_string(struct Stack*output_stack, struct Stack*local_stack,
    struct Number**factors, size_t factor_count)
{
    size_t char_count = 0;
    for (size_t i = 0; i < factor_count; ++i)
    {
        *(char*)ALLOCATE(output_stack, char) = '*';
        char_count += 1 + number_string(output_stack, local_stack, factors[i]);
    }
    return char_count;
}

size_t product_string(struct Stack*output_stack, struct Stack*local_stack, struct Number**factors,
    size_t factor_count)
{
    size_t char_count = number_string(output_stack, local_stack, factors[0]);
    return char_count +
        trailing_factor_string(output_stack, local_stack, factors + 1, factor_count - 1);
}

size_t number_string(struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        size_t char_count = integer_string(output_stack, local_stack, a->value->numerator);
        if (!integer_equals(a->value->denominator, &one))
        {
            *(char*)ALLOCATE(output_stack, char) = '/';
            char_count += 1 + integer_string(output_stack, local_stack, a->value->denominator);
        }
        return char_count;
    }
    case '^':
    {
        size_t char_count;
        if (a->radicand->operation != 'r' || a->radicand->value->numerator->sign < 0)
        {
            *(char*)ALLOCATE(output_stack, char) = '(';
            char_count = 2 + number_string(output_stack, local_stack, a->radicand);
            *(char*)ALLOCATE(output_stack, char) = ')';
        }
        else
        {
            char_count = number_string(output_stack, local_stack, a->radicand);
        }
        char*exponent = ARRAY_ALLOCATE(output_stack, 4, char);
        exponent[0] = '^';
        exponent[1] = '(';
        exponent[2] = '1';
        exponent[3] = '/';
        char_count += 5 + integer_string(output_stack, local_stack, a->index);
        *(char*)ALLOCATE(output_stack, char) = ')';
        return char_count;
    }
    case '*':
    {
        size_t char_count;
        if (a->elements[0]->operation == 'r')
        {
            if (integer_equals(a->elements[0]->value->numerator, &INT(1, -)))
            {
                *(char*)ALLOCATE(output_stack, char) = '-';
                char_count = 1 + product_string(output_stack, local_stack, a->elements + 1,
                    a->element_count - 1);
            }
            else if (integer_equals(a->elements[0]->value->numerator, &one))
            {
                char_count = product_string(output_stack, local_stack, a->elements + 1,
                    a->element_count - 1);
            }
            else
            {
                char_count =
                    integer_string(output_stack, local_stack, a->elements[0]->value->numerator);
                char*next_factor = local_stack->cursor;
                size_t next_factor_char_count =
                    number_string(local_stack, output_stack, a->elements[1]);
                if (next_factor[0] != '(')
                {
                    *(char*)ALLOCATE(output_stack, char) = '*';
                    char_count += 1;
                }
                char*next_factor_out = ARRAY_ALLOCATE(output_stack, next_factor_char_count, char);
                memcpy(next_factor_out, next_factor, next_factor_char_count);
                char_count += next_factor_char_count + trailing_factor_string(output_stack,
                    local_stack, a->elements + 2, a->element_count - 2);
                local_stack->cursor = next_factor;
            }
            if (!integer_equals(a->elements[0]->value->denominator, &one))
            {
                *(char*)ALLOCATE(output_stack, char) = '/';
                char_count += 1 + integer_string(output_stack, local_stack,
                    a->elements[0]->value->denominator);
            }
        }
        else
        {
            char_count = product_string(output_stack, local_stack, a->elements, a->element_count);
        }
        return char_count;
    }
    case '+':
    {
        size_t char_count = number_string(output_stack, local_stack, a->elements[0]);
        for (size_t i = 1; i < a->element_count; ++i)
        {
            char*sum_string = local_stack->cursor;
            size_t sum_char_count = number_string(local_stack, output_stack, a->elements[i]);
            if (sum_string[0] != '-')
            {
                *(char*)ALLOCATE(output_stack, char) = '+';
                char_count += 1;
            }
            char*sum_string_out = ARRAY_ALLOCATE(output_stack, sum_char_count, char);
            memcpy(sum_string_out, sum_string, sum_char_count);
            char_count += sum_char_count;
            local_stack->cursor = sum_string;
        }
        return char_count;
    }
    default:
        crash("Number operation not recognized.");
    }
}