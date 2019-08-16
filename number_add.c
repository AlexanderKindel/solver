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

void sun_convert_terms_to_new_generator(struct Stack*output_stack, struct Stack*local_stack,
    struct Matrix*out, struct Number*a, struct RationalPolynomial*old_generator_in_terms_of_new,
    struct RationalPolynomial*new_generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    out->width = a->element_count;
    out->height = 0;
    size_t max_term_coefficient_count = 0;
    for (size_t i = 0; i < a->element_count; ++i)
    {
        max_term_coefficient_count = max(max_term_coefficient_count,
            a->terms_in_terms_of_generator[i]->coefficient_count);
    }
    struct RationalPolynomial**old_generator_powers =
        ARRAY_ALLOCATE(local_stack, max_term_coefficient_count, struct RationalPolynomial*);
    old_generator_powers[0] = &rational_polynomial_one;
    for (size_t i = 1; i < max_term_coefficient_count; ++i)
    {
        old_generator_powers[i] = number_field_element_multiply(local_stack, output_stack,
            old_generator_powers[i - 1], old_generator_in_terms_of_new,
            new_generator_minimal_polynomial);
        out->height = max(out->height, old_generator_powers[i]->coefficient_count);
    }
    out->rows = ARRAY_ALLOCATE(output_stack, out->height, struct Rational**);
    for (size_t i = 0; i < out->height; ++i)
    {
        out->rows[i] = ARRAY_ALLOCATE(output_stack, out->width, struct Rational*);
    }
    for (size_t i = 0; i < out->width; ++i)
    {
        struct RationalPolynomial*term_in_terms_of_new_generator =
            (struct RationalPolynomial*)&polynomial_zero;
        for (size_t j = 0; j < a->terms_in_terms_of_generator[i]->coefficient_count; ++j)
        {
            term_in_terms_of_new_generator = rational_polynomial_add(local_stack, output_stack,
                rational_polynomial_rational_multiply(local_stack, output_stack,
                    old_generator_powers[j], a->terms_in_terms_of_generator[i]->coefficients[j]),
                term_in_terms_of_new_generator);
        }
        for (size_t j = 0; j < term_in_terms_of_new_generator->coefficient_count; ++j)
        {
            out->rows[j][i] =
                rational_copy(output_stack, term_in_terms_of_new_generator->coefficients[j]);
        }
        for (size_t j = term_in_terms_of_new_generator->coefficient_count; j < out->height; ++j)
        {
            out->rows[j][i] = &rational_zero;
        }
    }
}

void sum_append_term(struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Number*new_term, struct RationalPolynomial*new_term_in_terms_of_generator)
{
    a->elements[a->element_count - 1] = number_copy(output_stack, new_term);
    a->terms_in_terms_of_generator[a->element_count - 1] =
        rational_polynomial_copy(output_stack, new_term_in_terms_of_generator);
    for (size_t i = a->element_count; i-- > 1;)
    {
        if (number_formal_compare(output_stack, local_stack, a->elements[i],
            a->elements[i - 1]) <= 0)
        {
            return;
        }
        POINTER_SWAP(a->elements[i], a->elements[i - 1]);
        POINTER_SWAP(a->terms_in_terms_of_generator[i], a->terms_in_terms_of_generator[i - 1]);
    }
}

struct Number*sum_incorporate_term_in_terms_of_new_generator(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Number*new_term,
    struct Number*new_generator, struct Matrix*a_terms_in_terms_of_new_generator,
    struct RationalPolynomial*new_term_in_terms_of_new_generator)
{
    void*output_stack_savepoint = output_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '+';
    if (new_term_in_terms_of_new_generator->coefficient_count >
        a_terms_in_terms_of_new_generator->height)
    {
        out->generator = number_copy(output_stack, a->generator);
        out->element_count = a->element_count + 1;
        out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(output_stack, out->element_count, struct RationalPolynomial*);
        for (size_t i = 0; i < a->element_count; ++i)
        {
            out->elements[i] = number_copy(output_stack, a->elements[i]);
            size_t highest_nonzero_coefficient_index = a_terms_in_terms_of_new_generator->height;
            while (highest_nonzero_coefficient_index > 0)
            {
                --highest_nonzero_coefficient_index;
                if (a_terms_in_terms_of_new_generator->rows
                    [highest_nonzero_coefficient_index][i]->numerator->value_count)
                {
                    break;
                }
            }
            out->terms_in_terms_of_generator[i] =
                polynomial_allocate(output_stack, highest_nonzero_coefficient_index + 1);
            for (size_t j = 0; j <= highest_nonzero_coefficient_index; ++j)
            {
                out->terms_in_terms_of_generator[i]->coefficients[j] = rational_copy(output_stack,
                    a_terms_in_terms_of_new_generator->rows[highest_nonzero_coefficient_index][j]);
            }
        }
        sum_append_term(output_stack, local_stack, out, new_term,
            new_term_in_terms_of_new_generator);
    }
    else
    {
        void*local_stack_savepoint = output_stack->cursor;
        out->generator = number_copy(output_stack, new_generator);
        out->elements = ARRAY_ALLOCATE(output_stack, a->element_count + 1, struct Number*);
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(output_stack, a->element_count + 1, struct RationalPolynomial*);
        struct Rational**augmentation = ARRAY_ALLOCATE(output_stack,
            a_terms_in_terms_of_new_generator->height, struct Rational*);
        memcpy(augmentation, new_term_in_terms_of_new_generator->coefficients,
            new_term_in_terms_of_new_generator->coefficient_count * sizeof(struct Rational*));
        for (size_t i = new_term_in_terms_of_new_generator->coefficient_count;
            i < a_terms_in_terms_of_new_generator->height; ++i)
        {
            augmentation[i] = &rational_zero;
        }
        matrix_row_echelon_form(rational_multiply, rational_subtract, local_stack, output_stack,
            a_terms_in_terms_of_new_generator, augmentation);
        for (size_t i = a_terms_in_terms_of_new_generator->height;
            i-- > a_terms_in_terms_of_new_generator->width;)
        {
            if (augmentation[i]->numerator->value_count != 0)
            {
                memcpy(out->elements, a->elements, a->element_count * sizeof(struct Number*));
                out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
                    a->minimal_polynomial, new_term->minimal_polynomial);
                sum_append_term(output_stack, local_stack, out, new_term,
                    new_term_in_terms_of_new_generator);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
            else
            {
                size_t terms_remaining = a_terms_in_terms_of_new_generator->height - i - 1;
                memcpy(augmentation + i, augmentation + i + 1,
                    terms_remaining * sizeof(struct Rational*));
                memcpy(a_terms_in_terms_of_new_generator->rows + i,
                    a_terms_in_terms_of_new_generator->rows + i + 1,
                    terms_remaining * sizeof(struct Rational**));
                --a_terms_in_terms_of_new_generator->height;
            }
        }
        matrix_diagonalize(local_stack, output_stack, a_terms_in_terms_of_new_generator,
            augmentation);
        out->element_count = 0;
        for (size_t i = 0; i < a_terms_in_terms_of_new_generator->width; ++i)
        {
            struct Rational*scale =
                rational_add(local_stack, output_stack, augmentation[i], &rational_one);
            if (scale->numerator->value_count)
            {
                out->elements[out->element_count] =
                    number_rational_multiply(output_stack, local_stack, a->elements[i], scale);
                out->terms_in_terms_of_generator[out->element_count] =
                    rational_polynomial_rational_multiply(output_stack, local_stack,
                        a->terms_in_terms_of_generator[i], scale);
                ++out->element_count;
            }
        }
        switch (out->element_count)
        {
        case 0:
        {
            output_stack->cursor = output_stack_savepoint;
            local_stack->cursor = local_stack_savepoint;
            return number_rational_initialize(output_stack, &rational_zero);
        }
        case 1:
        {
            memmove(output_stack_savepoint, out->elements[0], sizeof(struct Number));
            output_stack->cursor = (void*)((size_t)output_stack_savepoint + sizeof(struct Number));
            local_stack->cursor = local_stack_savepoint;
            return output_stack_savepoint;
        }
        }
        local_stack->cursor = local_stack_savepoint;
    }
    out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
        a->minimal_polynomial, new_term->minimal_polynomial);
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
    struct NestedPolynomial*nested_minimal_polynomial =
        rational_polynomial_to_nested_polynomial(local_stack, a->minimal_polynomial);
    struct NestedPolynomial**minimal_polynomial_factors = ARRAY_ALLOCATE(local_stack,
        nested_minimal_polynomial->coefficient_count - 1, struct NestedPolynomial*);
    size_t candidate_factor_count = number_field_polynomial_factor(local_stack, output_stack,
        nested_minimal_polynomial, b->minimal_polynomial, minimal_polynomial_factors);
    struct RationalPolynomial**candidate_factors =
        (struct RationalPolynomial**)minimal_polynomial_factors;
    for (size_t i = 0; i < candidate_factor_count;)
    {
        if (minimal_polynomial_factors[i]->coefficient_count == 2)
        {
            candidate_factors[i] = rational_polynomial_negative(local_stack,
                minimal_polynomial_factors[i]->coefficients[0]);
            ++i;
        }
        else
        {
            --candidate_factor_count;
            minimal_polynomial_factors[i] = minimal_polynomial_factors[candidate_factor_count];
        }
    }
    struct Number**conjugates = number_conjugates(output_stack, local_stack, a);
    size_t conjugate_count = a->minimal_polynomial->coefficient_count - 2;
    conjugates[0] = conjugates[conjugate_count];
    struct RectangularEstimate*conjugate_estimates =
        ARRAY_ALLOCATE(local_stack, conjugate_count, struct RectangularEstimate*);
    for (size_t i = 0; i < conjugate_count; ++i)
    {
        number_rectangular_estimate(local_stack, output_stack, conjugate_estimates + i,
            conjugates[i], &rational_one);
    }
    size_t*conjugate_indices = ARRAY_ALLOCATE(local_stack, conjugate_count, size_t);
    struct RectangularEstimate a_estimate;
    number_rectangular_estimate(local_stack, output_stack, &a_estimate, a, &rational_one);
    struct RectangularEstimate b_estimate;
    number_rectangular_estimate(local_stack, output_stack, &b_estimate, b, &rational_one);
    for (size_t i = 0; i < candidate_factor_count; ++i)
    {
        size_t conjugate_index_count = conjugate_count;
        for (size_t j = 0; j < conjugate_count; ++j)
        {
            conjugate_indices[j] = j;
        }
        while (true)
        {
            struct RectangularEstimate factor_at_b_estimate;
            rational_polynomial_evaluate_at_rectangular_estimate(local_stack, output_stack,
                &factor_at_b_estimate, candidate_factors[i], &b_estimate);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &a_estimate,
                &factor_at_b_estimate))
            {
                break;
            }
            for (size_t j = 0; j < conjugate_index_count;)
            {
                struct RectangularEstimate*conjugate_estimate = conjugate_estimates + j;
                if (rectangular_estimates_are_disjoint(output_stack, local_stack,
                    conjugate_estimate, &factor_at_b_estimate))
                {
                    --conjugate_index_count;
                    if (conjugate_index_count == 0)
                    {
                        struct RationalPolynomial*out =
                            rational_polynomial_copy(output_stack, candidate_factors[i]);
                        local_stack->cursor = local_stack_savepoint;
                        return out;
                    }
                    conjugate_indices[j] = conjugate_indices[conjugate_index_count];
                }
                else
                {
                    ++j;
                }
                number_rectangular_estimate_halve_dimensions(local_stack, output_stack,
                    conjugate_estimate, conjugates[conjugate_indices[j]]);
            }
            number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &a_estimate,
                a);
            number_rectangular_estimate_halve_dimensions(local_stack, output_stack, &a_estimate,
                b);
        }
    }
    return 0;
}

struct Number*sum_incorporate_term(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Number*new_term)
{
    void*local_stack_savepoint = local_stack->cursor;
    for (size_t i = 0; i < a->element_count; ++i)
    {
        struct Rational*new_term_in_terms_of_old;
        switch (new_term->operation)
        {
        case 'r':
            if (a->elements[i]->operation == 'r')
            {
                new_term_in_terms_of_old = rational_divide(local_stack, output_stack,
                    new_term->value, a->elements[i]->value);
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
                    a->elements[i]) == 0)
                {
                    new_term_in_terms_of_old = new_term->elements[0]->value;
                    break;
                }
            }
        case '^':
            if (number_formal_compare(output_stack, local_stack, new_term, a->elements[i]) == 0)
            {
                new_term_in_terms_of_old = &rational_one;
                break;
            }
            if (a->elements[i]->operation == '*' && a->elements[i]->elements[0]->operation == 'r')
            {
                struct Number*right_factors = number_copy(local_stack, a->elements[i]->elements[0]);
                ++right_factors->elements;
                --right_factors->element_count;
                if (!number_formal_compare(output_stack, local_stack, new_term, right_factors))
                {
                    new_term_in_terms_of_old =
                        rational_reciprocal(local_stack, a->elements[i]->elements[0]->value);
                    break;
                }
            }
            continue;
        }
        struct Number*out;
        if (rational_equals(new_term_in_terms_of_old, &(struct Rational){&INT(1, -), &one}))
        {
            if (a->element_count == 2)
            {
                out = number_copy(output_stack, a->elements[a->element_count - i - 1]);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
            out = ALLOCATE(output_stack, struct Number);
            out->element_count = a->element_count - 1;
            out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, out->element_count, struct RationalPolynomial*);
        }
        else
        {
            out = ALLOCATE(output_stack, struct Number);
            out->element_count = a->element_count;
            out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, out->element_count, struct RationalPolynomial*);
            out->elements[i] = number_rational_multiply(output_stack, local_stack, a->elements[i],
                rational_integer_add(local_stack, output_stack, new_term_in_terms_of_old,
                    &one));
            out->terms_in_terms_of_generator[i] =
                rational_polynomial_rational_multiply(output_stack, local_stack,
                    a->terms_in_terms_of_generator[i], new_term_in_terms_of_old);
        }
        out->operation = '+';
        out->generator = number_copy(output_stack, a->generator);
        for (size_t j = 0; j < i; ++j)
        {
            out->elements[j] = number_copy(output_stack, a->elements[j]);
            out->terms_in_terms_of_generator[j] =
                rational_polynomial_copy(output_stack, a->terms_in_terms_of_generator[j]);
        }
        for (size_t j = 1; j < a->element_count - i; ++j)
        {
            size_t out_term_index = out->element_count - j;
            size_t a_term_index = a->element_count - j;
            out->elements[out_term_index] = number_copy(output_stack, a->elements[a_term_index]);
            out->terms_in_terms_of_generator[out_term_index] =
                rational_polynomial_copy(output_stack,
                    a->terms_in_terms_of_generator[a_term_index]);
        }
        out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
            a->minimal_polynomial, new_term->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct RationalPolynomial*new_term_in_terms_of_a_generator =
        number_a_in_terms_of_b(local_stack, output_stack, new_term, a->generator);
    if (new_term_in_terms_of_a_generator)
    {
        struct Matrix a_terms_in_terms_of_a_generator = { ARRAY_ALLOCATE(local_stack,
            a->generator->minimal_polynomial->coefficient_count - 1, struct Rational**),
            a->element_count, a->generator->minimal_polynomial->coefficient_count - 1 };
        for (size_t i = 0; i < a_terms_in_terms_of_a_generator.height; ++i)
        {
            a_terms_in_terms_of_a_generator.rows[i] = ARRAY_ALLOCATE(local_stack,
                a_terms_in_terms_of_a_generator.width, struct Rational*);
        }
        size_t max_term_coefficient_count = 0;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            max_term_coefficient_count = max(max_term_coefficient_count,
                a->terms_in_terms_of_generator[i]->coefficient_count);
            for (size_t j = 0; j < a->terms_in_terms_of_generator[i]->coefficient_count; ++j)
            {
                a_terms_in_terms_of_a_generator.rows[j][i] =
                    a->terms_in_terms_of_generator[i]->coefficients[j];
            }
            for (size_t j = a->terms_in_terms_of_generator[i]->coefficient_count;
                j < a_terms_in_terms_of_a_generator.height; ++j)
            {
                a_terms_in_terms_of_a_generator.rows[j][i] = &rational_zero;
            }
        }
        a_terms_in_terms_of_a_generator.height = max_term_coefficient_count;
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a, new_term, a->generator, &a_terms_in_terms_of_a_generator,
            new_term_in_terms_of_a_generator);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct RationalPolynomial x = { 2, &rational_zero, &rational_one };
    struct RationalPolynomial*old_generator_in_terms_of_new_term =
        number_a_in_terms_of_b(local_stack, output_stack, a->generator, new_term);
    if (old_generator_in_terms_of_new_term)
    {
        struct Matrix a_terms_in_terms_of_new_term;
        sun_convert_terms_to_new_generator(local_stack, output_stack, &a_terms_in_terms_of_new_term,
            a, old_generator_in_terms_of_new_term, new_term->minimal_polynomial);
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a, new_term, new_term, &a_terms_in_terms_of_new_term, &x);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Number*new_generator = ALLOCATE(local_stack, struct Number);
    new_generator->operation = 'g';
    new_generator->element_count = a->generator->element_count + 1;
    new_generator->elements =
        ARRAY_ALLOCATE(local_stack, new_generator->element_count, struct Number*);
    memcpy(new_generator->elements, a->generator->elements,
        a->generator->element_count * sizeof(struct Number*));
    struct Rational*k = &rational_one;
    while (true)
    {
        new_generator->elements[a->generator->element_count] = number_rational_multiply(local_stack,
            output_stack, new_term, k);
        new_generator->minimal_polynomial = sum_minimal_polynomial(local_stack, output_stack,
            new_generator, a->generator->minimal_polynomial,
            new_generator->elements[a->generator->element_count]->minimal_polynomial);
        struct RationalPolynomial*term_in_terms_of_new_generator =
            number_a_in_terms_of_b(local_stack, output_stack, new_term, new_generator);
        if (term_in_terms_of_new_generator)
        {
            struct Matrix a_terms_in_terms_of_new_generator;
            sun_convert_terms_to_new_generator(local_stack, output_stack,
                &a_terms_in_terms_of_new_generator, a,
                rational_polynomial_subtract(local_stack, output_stack,
                    (struct RationalPolynomial*)&x,
                        rational_polynomial_rational_multiply(local_stack, output_stack,
                            term_in_terms_of_new_generator, k)), new_generator->minimal_polynomial);
            struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
                local_stack, a, new_term, new_generator, &a_terms_in_terms_of_new_generator,
                term_in_terms_of_new_generator);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        k->numerator = integer_add(local_stack, k->numerator, &one);
    }
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
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*a_term = ALLOCATE(local_stack, struct Number);
            a_term->operation = '+';
            a_term->generator = number_copy(local_stack, a);
            a_term->element_count = 1;
            a_term->elements = ALLOCATE(local_stack, struct Number*);
            a_term->elements[0] = a;
            a_term->terms_in_terms_of_generator =
                ALLOCATE(local_stack, struct RationalPolynomial*);
            a_term->terms_in_terms_of_generator[0] =
                &(struct RationalPolynomial){ 2, &rational_zero, &rational_one };
            a_term->minimal_polynomial = a->minimal_polynomial;
            struct Number*out = sum_incorporate_term(output_stack, local_stack, a_term, b);
            local_stack->cursor = local_stack_savepoint;
            return out;
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
            return sum_incorporate_term(output_stack, local_stack, a, b);
        }
    }
    return number_add(output_stack, local_stack, b, a);
}