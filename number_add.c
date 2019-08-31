#include "declarations.h"

struct RationalPolynomial*sum_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*power = &nested_polynomial_one;
    struct NestedPolynomial d = { 2, &(struct RationalPolynomial){2, &rational_zero, &rational_one},
        &(struct RationalPolynomial){1, &(struct Rational){&INT(1, -), &one}} };
    struct NestedPolynomial*t = (struct NestedPolynomial*)&polynomial_zero;
    for (size_t i = 0; i < left_minimal_polynomial->coefficient_count; ++i)
    {
        t = nested_polynomial_add(local_stack, output_stack, t,
            nested_polynomial_rational_polynomial_multiply(local_stack, output_stack, power,
                &(struct RationalPolynomial){ 1, left_minimal_polynomial->coefficients[i] }));
        power = nested_polynomial_multiply(local_stack, output_stack, power, &d);
    }
    struct RationalPolynomial*out =
        number_minimal_polynomial_from_annulling_polynomial(output_stack, local_stack,
            nested_polynomial_resultant(local_stack, output_stack, t,
                rational_polynomial_to_nested_polynomial(local_stack, right_minimal_polynomial)),
            a);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial**sum_convert_terms_to_new_generator(struct Stack*output_stack,
    struct Stack*local_stack, size_t term_count,
    struct RationalPolynomial**terms_in_terms_of_old_generator,
    struct RationalPolynomial*old_generator_in_terms_of_new,
    struct RationalPolynomial*new_generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial**out =
        ARRAY_ALLOCATE(output_stack, term_count, struct RationalPolynomial*);
    size_t max_term_coefficient_count = 0;
    for (size_t i = 0; i < term_count; ++i)
    {
        max_term_coefficient_count = max(max_term_coefficient_count,
            terms_in_terms_of_old_generator[i]->coefficient_count);
    }
    struct RationalPolynomial**old_generator_powers =
        ARRAY_ALLOCATE(local_stack, max_term_coefficient_count, struct RationalPolynomial*);
    old_generator_powers[0] = &rational_polynomial_one;
    for (size_t i = 1; i < max_term_coefficient_count; ++i)
    {
        old_generator_powers[i] = number_field_element_multiply(local_stack, output_stack,
            old_generator_powers[i - 1], old_generator_in_terms_of_new,
            new_generator_minimal_polynomial);
    }
    for (size_t i = 0; i < term_count; ++i)
    {
        out[i] = (struct RationalPolynomial*)&polynomial_zero;
        for (size_t j = 0; j < terms_in_terms_of_old_generator[i]->coefficient_count; ++j)
        {
            out[i] = rational_polynomial_add(local_stack, output_stack,
                rational_polynomial_rational_multiply(local_stack, output_stack,
                    old_generator_powers[j], terms_in_terms_of_old_generator[i]->coefficients[j]),
                out[i]);
        }
        out[i] = rational_polynomial_copy(output_stack, out[i]);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_a_in_terms_of_b(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Number*b)
{
    if ((b->minimal_polynomial->coefficient_count - 1) %
        (a->minimal_polynomial->coefficient_count - 1) != 0)
    {
        return 0;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*nested_a_minimal_polynomial =
        rational_polynomial_to_nested_polynomial(local_stack,
            rational_polynomial_euclidean_remainder(local_stack, output_stack,
                a->minimal_polynomial, b->minimal_polynomial));
    switch (nested_a_minimal_polynomial->coefficient_count)
    {
    case 0:
        struct RationalPolynomial*out = polynomial_allocate(output_stack, 2);
        out->coefficients[0] = &rational_zero;
        out->coefficients[1] = &rational_one;
        return out;
    case 1:
        local_stack->cursor = local_stack_savepoint;
        return 0;
    }
    struct NestedPolynomial**a_minimal_polynomial_factors = ARRAY_ALLOCATE(local_stack,
        nested_a_minimal_polynomial->coefficient_count - 1, struct NestedPolynomial*);
    size_t candidate_factor_count = number_field_polynomial_factor(local_stack, output_stack,
        nested_a_minimal_polynomial, b->minimal_polynomial, a_minimal_polynomial_factors);
    struct RationalPolynomial**candidate_factors =
        (struct RationalPolynomial**)a_minimal_polynomial_factors;
    for (size_t i = 0; i < candidate_factor_count;)
    {
        if (a_minimal_polynomial_factors[i]->coefficient_count == 2)
        {
            candidate_factors[i] = rational_polynomial_negative(local_stack,
                a_minimal_polynomial_factors[i]->coefficients[0]);
            ++i;
        }
        else
        {
            --candidate_factor_count;
            a_minimal_polynomial_factors[i] = a_minimal_polynomial_factors[candidate_factor_count];
        }
    }
    struct RectangularEstimate a_estimate;
    number_rectangular_estimate(local_stack, output_stack, &a_estimate, a, &rational_one);
    while (rational_polynomial_root_count_in_rectangle(local_stack, output_stack,
        a->minimal_polynomial,
        float_interval_to_rational_interval(local_stack, output_stack,
            a_estimate.real_part_estimate),
        float_interval_to_rational_interval(local_stack, output_stack,
            a_estimate.imaginary_part_estimate)) > 1)
    {
        number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &a_estimate, a);
    }
    struct RectangularEstimate b_estimate;
    number_rectangular_estimate(local_stack, output_stack, &b_estimate, b, &rational_one);
    for (size_t i = 0; i < candidate_factor_count; ++i)
    {
        struct Rational interval_size_for_evaluation =
        { &one, integer_from_size_t(local_stack, candidate_factors[i]->coefficient_count) };
        struct RectangularEstimate factor_at_b;
        while (true)
        {
            rational_polynomial_evaluate_at_rectangular_estimate(local_stack, output_stack,
                &factor_at_b, candidate_factors[i], &b_estimate, &interval_size_for_evaluation);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &a_estimate,
                &factor_at_b))
            {
                goto candidate_rejected;
            }
            if (rational_polynomial_root_count_in_rectangle(local_stack, output_stack,
                a->minimal_polynomial,
                float_interval_to_rational_interval(local_stack, output_stack,
                    factor_at_b.real_part_estimate),
                float_interval_to_rational_interval(local_stack, output_stack,
                    factor_at_b.imaginary_part_estimate)) == 1)
            {
                if (rectangular_estimate_a_contains_b(output_stack, local_stack, &a_estimate,
                    &factor_at_b))
                {
                    struct RationalPolynomial*out =
                        rational_polynomial_copy(output_stack, candidate_factors[i]);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
                break;
            }
            number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &b_estimate, b);
            interval_size_for_evaluation.denominator =
                integer_doubled(local_stack, interval_size_for_evaluation.denominator);
        }
        while (true)
        {
            if (rectangular_estimate_a_contains_b(output_stack, local_stack, &factor_at_b,
                &a_estimate))
            {
                struct RationalPolynomial*out =
                    rational_polynomial_copy(output_stack, candidate_factors[i]);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
            number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &a_estimate, a);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &a_estimate,
                &factor_at_b))
            {
                goto candidate_rejected;
            }
        }
    candidate_rejected:;
    }
    local_stack->cursor = local_stack_savepoint;
    return 0;
}

struct Number*sum_append_term(struct Stack*output_stack, struct Stack*local_stack,
    struct Number**old_terms, size_t old_term_count,
    struct RationalPolynomial**old_terms_in_terms_of_new_generator, struct Number*new_term,
    struct Number*new_generator, struct RationalPolynomial*new_term_in_terms_of_new_generator)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '+';
    out->element_count = old_term_count + 1;
    out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
    out->generator = number_copy(output_stack, new_generator);
    out->terms_in_terms_of_generator =
        ARRAY_ALLOCATE(output_stack, out->element_count, struct RationalPolynomial*);
    for (size_t i = 0; i < old_term_count; ++i)
    {
        out->elements[i] = number_copy(output_stack, old_terms[i]);
        out->terms_in_terms_of_generator[i] =
            rational_polynomial_copy(output_stack, old_terms_in_terms_of_new_generator[i]);
    }
    out->elements[old_term_count] = number_copy(output_stack, new_term);
    out->terms_in_terms_of_generator[old_term_count] =
        rational_polynomial_copy(output_stack, new_term_in_terms_of_new_generator);
    for (size_t i = out->element_count; i-- > 1;)
    {
        if (number_formal_compare(output_stack, local_stack, out->elements[i],
            out->elements[i - 1]) <= 0)
        {
            return out;
        }
        POINTER_SWAP(out->elements[i], out->elements[i - 1]);
        POINTER_SWAP(out->terms_in_terms_of_generator[i], out->terms_in_terms_of_generator[i - 1]);
    }
    return out;
}

struct Number*sum_incorporate_term_in_terms_of_new_generator(struct Stack*output_stack,
    struct Stack*local_stack, struct Number**a_terms, size_t a_term_count,
    struct RationalPolynomial*a_minimal_polynomial, struct Number*new_term,
    struct Number*new_generator, struct RationalPolynomial**a_terms_in_terms_of_new_generator,
    struct RationalPolynomial*new_term_in_terms_of_new_generator)
{
    void*local_stack_savepoint = output_stack->cursor;
    size_t max_term_coefficient_count = 0;
    for (size_t i = 0; i < a_term_count; ++i)
    {
        max_term_coefficient_count = max(max_term_coefficient_count,
            a_terms_in_terms_of_new_generator[i]->coefficient_count);
    }
    struct Matrix matrix =
    { ARRAY_ALLOCATE(local_stack, max_term_coefficient_count, struct Rational**),
        a_term_count, max_term_coefficient_count };
    for (size_t i = 0; i < matrix.height; ++i)
    {
        matrix.rows[i] = ARRAY_ALLOCATE(local_stack, matrix.width, struct Rational*);
    }
    for (size_t i = 0; i < a_term_count; ++i)
    {
        for (size_t j = 0; j < a_terms_in_terms_of_new_generator[i]->coefficient_count; ++j)
        {
            matrix.rows[j][i] = a_terms_in_terms_of_new_generator[i]->coefficients[j];
        }
        for (size_t j = a_terms_in_terms_of_new_generator[i]->coefficient_count; j < matrix.height;
            ++j)
        {
            matrix.rows[j][i] = &rational_zero;
        }
    }
    if (new_term_in_terms_of_new_generator->coefficient_count > matrix.height)
    {
        struct Number*out = sum_append_term(output_stack, local_stack, a_terms, a_term_count,
            a_terms_in_terms_of_new_generator, new_term, new_generator,
            new_term_in_terms_of_new_generator);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Rational**augmentation = ARRAY_ALLOCATE(local_stack, matrix.height, struct Rational*);
    memcpy(augmentation, new_term_in_terms_of_new_generator->coefficients,
        new_term_in_terms_of_new_generator->coefficient_count * sizeof(struct Rational*));
    for (size_t i = new_term_in_terms_of_new_generator->coefficient_count; i < matrix.height; ++i)
    {
        augmentation[i] = &rational_zero;
    }
    matrix_row_echelon_form(local_stack, output_stack, &matrix, augmentation);
    for (size_t i = matrix.height; i-- > matrix.width;)
    {
        if (augmentation[i]->numerator->value_count != 0)
        {
            struct Number*out = sum_append_term(output_stack, local_stack, a_terms, a_term_count,
                a_terms_in_terms_of_new_generator, new_term, new_generator,
                new_term_in_terms_of_new_generator);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        else
        {
            size_t terms_remaining = matrix.height - i - 1;
            memmove(augmentation + i, augmentation + i + 1,
                terms_remaining * sizeof(struct Rational*));
            memmove(matrix.rows + i, matrix.rows + i + 1,
                terms_remaining * sizeof(struct Rational**));
            --matrix.height;
        }
    }
    matrix_diagonalize(local_stack, output_stack, &matrix, augmentation);
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '+';
    out->generator = number_copy(output_stack, new_generator);
    out->elements = ARRAY_ALLOCATE(output_stack, a_term_count + 1, struct Number*);
    out->terms_in_terms_of_generator =
        ARRAY_ALLOCATE(output_stack, a_term_count + 1, struct RationalPolynomial*);
    out->element_count = 0;
    for (size_t i = 0; i < a_term_count; ++i)
    {
        struct Rational*scale =
            rational_add(local_stack, output_stack, augmentation[i], &rational_one);
        if (scale->numerator->value_count)
        {
            out->elements[out->element_count] =
                number_rational_multiply(output_stack, local_stack, a_terms[i], scale);
            out->terms_in_terms_of_generator[out->element_count] =
                rational_polynomial_rational_multiply(output_stack, local_stack,
                    a_terms_in_terms_of_new_generator[i], scale);
            ++out->element_count;
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Can handle the degenerate case of a_term_count == 1, and returns a sum even in degenerate cases of
//fewer than two terms. Leaves it to the caller to fill out the minimal_polynomial field since it
//will need to be calculated only in nondegenerate cases.
struct Number*sum_incorporate_term(struct Stack*output_stack, struct Stack*local_stack,
    struct Number**a_terms, size_t a_term_count,
    struct RationalPolynomial**a_terms_in_terms_of_generator, struct Number*a_generator,
    struct RationalPolynomial*a_minimal_polynomial, struct Number*new_term)
{
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a_term_count; ++i)
    {
        struct Rational*new_term_in_terms_of_old;
        switch (new_term->operation)
        {
        case 'r':
            if (a_terms[i]->operation == 'r')
            {
                new_term_in_terms_of_old = rational_divide(local_stack, output_stack,
                    new_term->value, a_terms[i]->value);
                break;
            }
            continue;
        case '*':
            if (new_term->elements[0]->operation == 'r')
            {
                struct Number*right_factors = number_copy(local_stack, new_term);
                ++right_factors->elements;
                --right_factors->element_count;
                if (number_formal_compare(output_stack, local_stack, right_factors,
                    a_terms[i]) == 0)
                {
                    new_term_in_terms_of_old = new_term->elements[0]->value;
                    break;
                }
            }
        case '^':
            if (number_formal_compare(output_stack, local_stack, new_term, a_terms[i]) == 0)
            {
                new_term_in_terms_of_old = &rational_one;
                break;
            }
            if (a_terms[i]->operation == '*' && a_terms[i]->elements[0]->operation == 'r')
            {
                struct Number right_factors;
                memcpy(&right_factors, a_terms[i]->elements[0], sizeof(struct Number));
                ++right_factors.elements;
                --right_factors.element_count;
                if (!number_formal_compare(output_stack, local_stack, new_term, &right_factors))
                {
                    new_term_in_terms_of_old =
                        rational_reciprocal(local_stack, a_terms[i]->elements[0]->value);
                    break;
                }
            }
            continue;
        }
        struct Number*out = ALLOCATE(output_stack, struct Number);
        out->operation = '+';
        out->generator = number_copy(output_stack, a_generator);
        if (rational_equals(new_term_in_terms_of_old, &(struct Rational){&INT(1, -), &one}))
        {
            out->element_count = a_term_count - 1;
            out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, out->element_count, struct RationalPolynomial*);
        }
        else
        {
            out->element_count = a_term_count;
            out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, out->element_count, struct RationalPolynomial*);
            out->elements[i] = number_rational_multiply(output_stack, local_stack, a_terms[i],
                rational_integer_add(local_stack, output_stack, new_term_in_terms_of_old, &one));
            out->terms_in_terms_of_generator[i] =
                rational_polynomial_rational_multiply(output_stack, local_stack,
                    a_terms_in_terms_of_generator[i], new_term_in_terms_of_old);
        }
        for (size_t j = 0; j < i; ++j)
        {
            out->elements[j] = number_copy(output_stack, a_terms[j]);
            out->terms_in_terms_of_generator[j] =
                rational_polynomial_copy(output_stack, a_terms_in_terms_of_generator[j]);
        }
        for (size_t j = 1; j < a_term_count - i; ++j)
        {
            size_t out_term_index = out->element_count - j;
            size_t a_term_index = a_term_count - j;
            out->elements[out_term_index] = number_copy(output_stack, a_terms[a_term_index]);
            out->terms_in_terms_of_generator[out_term_index] =
                rational_polynomial_copy(output_stack,
                    a_terms_in_terms_of_generator[a_term_index]);
        }
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct RationalPolynomial*new_term_in_terms_of_a_generator =
        number_a_in_terms_of_b(local_stack, output_stack, new_term, a_generator);
    if (new_term_in_terms_of_a_generator)
    {
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a_terms, a_term_count, a_minimal_polynomial, new_term, a_generator,
            a_terms_in_terms_of_generator, new_term_in_terms_of_a_generator);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct RationalPolynomial x = { 2, &rational_zero, &rational_one };
    struct RationalPolynomial*old_generator_in_terms_of_new_term =
        number_a_in_terms_of_b(local_stack, output_stack, a_generator, new_term);
    if (old_generator_in_terms_of_new_term)
    {
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a_terms, a_term_count, a_minimal_polynomial, new_term, new_term,
            sum_convert_terms_to_new_generator(local_stack, output_stack, a_term_count,
                a_terms_in_terms_of_generator, old_generator_in_terms_of_new_term,
                new_term->minimal_polynomial), &x);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Number*new_generator = ALLOCATE(local_stack, struct Number);
    new_generator->operation = '+';
    new_generator->generator = 0;
    switch (a_generator->operation)
    {
    case '+':
        new_generator->element_count = a_generator->element_count + 1;
        new_generator->elements =
            ARRAY_ALLOCATE(local_stack, new_generator->element_count, struct Number*);
        memcpy(new_generator->elements, a_generator->elements,
            a_generator->element_count * sizeof(struct Number*));
        break;
    default:
        new_generator->element_count = 2;
        new_generator->elements = ARRAY_ALLOCATE(local_stack, 2, struct Number*);
        new_generator->elements[0] = a_generator;
    }
    struct Rational*k = &(struct Rational) { &one, &one };
    while (true)
    {
        void*local_stack_loop_savepoint = local_stack->cursor;
        new_generator->elements[new_generator->element_count - 1] =
            number_rational_multiply(local_stack, output_stack, new_term, k);
        new_generator->minimal_polynomial = sum_minimal_polynomial(local_stack, output_stack,
            new_generator, a_generator->minimal_polynomial,
            new_generator->elements[new_generator->element_count - 1]->minimal_polynomial);
        struct RationalPolynomial*term_in_terms_of_new_generator =
            number_a_in_terms_of_b(local_stack, output_stack, new_term, new_generator);
        if (term_in_terms_of_new_generator)
        {
            struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
                local_stack, a_terms, a_term_count, a_minimal_polynomial, new_term, new_generator,
                sum_convert_terms_to_new_generator(local_stack, output_stack, a_term_count,
                    a_terms_in_terms_of_generator,
                    rational_polynomial_subtract(local_stack, output_stack, &x,
                        rational_polynomial_rational_multiply(local_stack, output_stack,
                            term_in_terms_of_new_generator, k)),
                    new_generator->minimal_polynomial), term_in_terms_of_new_generator);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        local_stack->cursor = local_stack_loop_savepoint;
        k->numerator = integer_add(local_stack, k->numerator, &one);
    }
}

struct Number*number_incorporate_term(struct Stack*output_stack, struct Stack*local_stack,
    struct Number**a_terms, size_t a_term_count,
    struct RationalPolynomial**a_terms_in_terms_of_generator, struct Number*a_generator,
    struct RationalPolynomial*a_minimal_polynomial, struct Number*new_term)
{
    void*output_stack_savepoint = output_stack->cursor;
    struct Number*out = sum_incorporate_term(output_stack, local_stack, a_terms, a_term_count,
        a_terms_in_terms_of_generator, a_generator, a_minimal_polynomial, new_term);
    switch (out->element_count)
    {
    case 0:
    {
        output_stack->cursor = output_stack_savepoint;
        return number_rational_initialize(output_stack, &rational_zero);
    }
    case 1:
    {
        memmove(output_stack_savepoint, out->elements[0], sizeof(struct Number));
        output_stack->cursor = (void*)((size_t)output_stack_savepoint + sizeof(struct Number));
        return output_stack_savepoint;
    }
    }
    out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
        a_minimal_polynomial, new_term->minimal_polynomial);
    return out;
}

struct Number*number_add(struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value->numerator->value_count == 0)
        {
            return number_copy(output_stack, b);
        }
        if (b->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*out = number_rational_initialize(output_stack,
                rational_add(local_stack, output_stack, a->value, b->value));
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        break;
    case '^':
    case '*':
        switch (b->operation)
        {
        case 'r':
            struct Number*out = ALLOCATE(output_stack, struct Number);
            out->operation = '+';
            out->element_count = 2;
            out->generator = number_copy(output_stack, a);
            out->elements = ARRAY_ALLOCATE(output_stack, 2, struct Number);
            out->elements[0] = out->generator;
            out->elements[1] = number_copy(output_stack, b);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, 2, struct RationalPolynomial);
            out->terms_in_terms_of_generator[0] = polynomial_allocate(output_stack, 2);
            out->terms_in_terms_of_generator[0]->coefficients[0] = &rational_zero;
            out->terms_in_terms_of_generator[0]->coefficients[1] = &rational_one;
            out->terms_in_terms_of_generator[1] = polynomial_allocate(output_stack, 1);
            out->terms_in_terms_of_generator[1]->coefficients[0] =
                rational_copy(output_stack, b->value);
            out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
                a->minimal_polynomial, b->minimal_polynomial);
            return out;
        case '^':
        case '*':   
        {
            struct RationalPolynomial*x =
                &(struct RationalPolynomial) { 2, &rational_zero, &rational_one };
            return number_incorporate_term(output_stack, local_stack, &a, 1, &x, a,
                a->minimal_polynomial, b);
        }
        }
        break;
    case '+':
        if (b->operation == '+')
        {
            void*local_stack_savepoint = local_stack->cursor;
            for (size_t i = b->element_count; i-- > 1;)
            {
                a = number_add(local_stack, output_stack, a, b->elements[i]);
            }
            a = number_add(output_stack, local_stack, a, b->elements[0]);
            local_stack->cursor = local_stack_savepoint;
            return a;
        }
        else
        {
            return number_incorporate_term(output_stack, local_stack, a->elements, a->element_count,
                a->terms_in_terms_of_generator, a->generator, a->minimal_polynomial, b);
        }
    }
    return number_add(output_stack, local_stack, b, a);
}