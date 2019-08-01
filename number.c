#include "declarations.h"

struct Number*number_copy(struct Stack*output_stack, struct Number*a)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = a->operation;
    out->minimal_polynomial = rational_polynomial_copy(output_stack, a->minimal_polynomial);
    switch (a->operation)
    {
    case 'r':
        out->value.numerator = integer_copy(output_stack, a->value.numerator);
        out->value.denominator = integer_copy(output_stack, a->value.denominator);
        return out;
    case '^':
    case '*':
        out->left = number_copy(output_stack, a->left);
        out->right = number_copy(output_stack, a->right);
        return out;
    case '+':
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(output_stack, a->term_count, struct RationalPolynomial*);
        for (size_t i = 0; i < a->term_count; ++i)
        {
            out->terms_in_terms_of_generator[i] =
                rational_polynomial_copy(output_stack, a->terms_in_terms_of_generator[i]);
        }
    }
    out->term_count = a->term_count;
    out->terms = ARRAY_ALLOCATE(output_stack, out->term_count, struct Number*);
    for (size_t i = 0; i < a->term_count; ++i)
    {
        out->terms[i] = number_copy(output_stack, a->terms[i]);
    }
    return out;
}

struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = 'r';
    out->value.numerator = integer_copy(output_stack, value->numerator);
    out->value.denominator = integer_copy(output_stack, value->denominator);
    out->minimal_polynomial = rational_minimal_polynomial(output_stack, value);
    return out;
}

//Assigns radicand and exponent directly to the left and right fields of the return value without
//copying them to output_stack.
struct Number*number_surd_initialize(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*radicand, struct Number*exponent)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '^';
    out->left = radicand;
    out->right = exponent;
    size_t surd_index = integer_to_size_t(exponent->value.denominator);
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
    struct RectangularEstimate zero_rectangular_estimate = { &zero_estimate, &zero_estimate };
    struct Rational*interval_size = &rational_one;
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct RectangularEstimate evaluation_estimate;
            rational_polynomial_estimate_evaluation(local_stack, output_stack, &evaluation_estimate,
                candidates[i], a, interval_size);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &evaluation_estimate,
                &zero_rectangular_estimate))
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
        interval_size->denominator = integer_doubled(local_stack, interval_size->denominator);
    }
}

struct RationalPolynomial*sum_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*power = &nested_polynomial_one;
    struct NestedPolynomial*d = polynomial_allocate(local_stack, 2);
    d->coefficients[0] = polynomial_allocate(local_stack, 2);
    d->coefficients[0]->coefficients[0] = &rational_zero;
    d->coefficients[0]->coefficients[1] = &rational_one;
    d->coefficients[1] = polynomial_allocate(local_stack, 1);
    d->coefficients[1]->coefficients[0]->numerator = integer_initialize(local_stack, 1, -1);
    d->coefficients[1]->coefficients[0]->denominator = &one;
    struct NestedPolynomial*t = polynomial_allocate(local_stack, 0);
    for (size_t i = 0; i < left_minimal_polynomial->coefficient_count; ++i)
    {
        struct RationalPolynomial*coefficient = polynomial_allocate(local_stack, 1);
        coefficient->coefficients[0] = left_minimal_polynomial->coefficients[i];
        t = nested_polynomial_add(local_stack, output_stack, t,
            nested_polynomial_rational_polynomial_multiply(local_stack, output_stack, power,
                coefficient));
        power = nested_polynomial_multiply(local_stack, output_stack, power, d);
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
    out->width = a->term_count;
    out->height = 0;
    size_t max_term_coefficient_count = 0;
    for (size_t i = 0; i < a->term_count; ++i)
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

struct Number*sum_incorporate_term_in_terms_of_new_generator(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Number*new_term,
    struct Number*new_generator, struct Matrix*a_terms_in_terms_of_new_generator,
    struct RationalPolynomial*new_term_in_terms_of_new_generator)
{
    struct Number*out;
    if (new_term_in_terms_of_new_generator->coefficient_count >
        a_terms_in_terms_of_new_generator->height)
    {
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '+';
        out->generator = number_copy(output_stack, a->generator);
        out->term_count = a->term_count + 1;
        out->terms = ARRAY_ALLOCATE(output_stack, out->term_count, struct Number*);
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(output_stack, out->term_count, struct RationalPolynomial*);
        for (size_t i = 0; i < a->term_count; ++i)
        {
            out->terms[i] = number_copy(output_stack, a->terms[i]);
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
        out->terms[a->term_count] = number_copy(output_stack, new_term);
        out->terms_in_terms_of_generator[a->term_count] =
            rational_polynomial_copy(output_stack, new_term_in_terms_of_new_generator);
        return out;
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        out = ALLOCATE(local_stack, struct Number);
        out->operation = '+';
        out->generator = number_copy(output_stack, new_generator);
        out->terms = ARRAY_ALLOCATE(local_stack, a->term_count + 1, struct Number*);
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(local_stack, a->term_count + 1, struct RationalPolynomial*);
        struct Rational**augmentation =
            ARRAY_ALLOCATE(local_stack, a_terms_in_terms_of_new_generator->height, struct Rational*);
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
                memcpy(out->terms, a->terms, a->term_count * sizeof(struct Number*));
                out->terms[a->term_count] = new_term;
                out->terms_in_terms_of_generator[a->term_count] =
                    new_term_in_terms_of_new_generator;
                out = number_copy(output_stack, out);
                out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
                    a->minimal_polynomial, new_term->minimal_polynomial);
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
        matrix_diagonalize(local_stack, output_stack, a_terms_in_terms_of_new_generator, augmentation);
        out->term_count = 0;
        for (size_t i = 0; i < a_terms_in_terms_of_new_generator->width; ++i)
        {
            struct Rational*scale =
                rational_add(local_stack, output_stack, augmentation[i], &rational_one);
            if (scale->numerator->value_count)
            {
                out->terms[out->term_count] =
                    number_rational_multiply(local_stack, output_stack, a->terms[i], scale);
                out->terms_in_terms_of_generator[out->term_count] =
                    rational_polynomial_rational_multiply(local_stack, output_stack,
                        a->terms_in_terms_of_generator[i], scale);
                ++out->term_count;
            }
        }
        switch (out->term_count)
        {
        case 0:
        {
            out = number_rational_initialize(output_stack, &rational_zero);
            break;
        }
        case 1:
        {
            out = number_copy(output_stack, out->terms[0]);
            break;
        }
        default:
            out = number_copy(output_stack, out);
            out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
                a->minimal_polynomial, new_term->minimal_polynomial);
        }
        local_stack->cursor = local_stack_savepoint;
    }
    return out;
}

struct RationalPolynomial*number_a_in_terms_of_b(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Number*b)
{
    if (a->operation == 'r')
    {
        struct RationalPolynomial*out = polynomial_allocate(output_stack, 1);
        out->coefficients[0] = rational_copy(output_stack, &a->value);
        return out;
    }
    if ((b->minimal_polynomial->coefficient_count - 1) %
        (a->minimal_polynomial->coefficient_count - 1) != 0)
    {
        return 0;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*nested_minimal_polynomial =
        rational_polynomial_to_nested_polynomial(local_stack, a->minimal_polynomial);
    struct NestedPolynomial**minimal_polynomial_factors = stack_slot_allocate(local_stack,
        (nested_minimal_polynomial->coefficient_count - 1) * sizeof(struct NestedPolynomial*),
        _Alignof(struct NestedPolynomial*));
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
    size_t*conjugate_indices = ARRAY_ALLOCATE(local_stack, conjugate_count, size_t);
    struct Rational*interval_size = &rational_one;
    for (size_t i = 0; i < candidate_factor_count; ++i)
    {
        size_t conjugate_index_count = conjugate_count;
        for (size_t j = 0; j < conjugate_count; ++j)
        {
            conjugate_indices[j] = j;
        }
        while (true)
        {
            struct RectangularEstimate a_estimate;
            number_rectangular_estimate(local_stack, output_stack, &a_estimate, a, interval_size);
            struct RectangularEstimate factor_at_b_estimate;
            rational_polynomial_estimate_evaluation(local_stack, output_stack,
                &factor_at_b_estimate, candidate_factors[i], b, interval_size);
            if (rectangular_estimates_are_disjoint(output_stack, local_stack, &a_estimate,
                &factor_at_b_estimate))
            {
                break;
            }
            for (size_t j = 0; j < conjugate_index_count;)
            {
                struct Number*conjugate = conjugates[conjugate_indices[j]];
                struct RectangularEstimate candidate_estimate;
                number_rectangular_estimate(local_stack, output_stack, &candidate_estimate,
                    conjugate, interval_size);
                if (rectangular_estimates_are_disjoint(output_stack, local_stack,
                    &candidate_estimate, &factor_at_b_estimate))
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
            }
            interval_size->denominator = integer_doubled(local_stack, interval_size->denominator);
        }
    }
    return 0;
}

struct Number*sum_incorporate_term(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Number*new_term)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*new_term_in_terms_of_a_generator =
        number_a_in_terms_of_b(local_stack, output_stack, new_term, a->generator);
    if (new_term_in_terms_of_a_generator)
    {
        struct Matrix a_terms_in_terms_of_a_generator = { ARRAY_ALLOCATE(local_stack,
            a->generator->minimal_polynomial->coefficient_count - 1, struct Rational**),
            a->term_count, a->generator->minimal_polynomial->coefficient_count - 1 };
        for (size_t i = 0; i < a_terms_in_terms_of_a_generator.height; ++i)
        {
            a_terms_in_terms_of_a_generator.rows[i] = ARRAY_ALLOCATE(local_stack, 
                a_terms_in_terms_of_a_generator.width, struct Rational*);
        }
        size_t max_term_coefficient_count = 0;
        for (size_t i = 0; i < a->term_count; ++i)
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
    struct { size_t coefficient_count; struct Rational*coefficients[2]; }x =
        { 2, { &rational_zero, &rational_one } };
    struct RationalPolynomial*old_generator_in_terms_of_new_term =
        number_a_in_terms_of_b(local_stack, output_stack, a->generator, new_term);
    if (old_generator_in_terms_of_new_term)
    {
        struct Matrix a_terms_in_terms_of_new_term;
        sun_convert_terms_to_new_generator(local_stack, output_stack, &a_terms_in_terms_of_new_term,
            a, old_generator_in_terms_of_new_term, new_term->minimal_polynomial);
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a, new_term, new_term, &a_terms_in_terms_of_new_term,
            (struct RationalPolynomial*)&x);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Number*new_generator = ALLOCATE(local_stack, struct Number);
    new_generator->operation = 'g';
    new_generator->term_count = a->generator->term_count + 1;
    new_generator->terms = ARRAY_ALLOCATE(local_stack, new_generator->term_count, struct Number*);
    memcpy(new_generator->terms, a->generator->terms,
        a->generator->term_count * sizeof(struct Number*));
    struct Rational*k = &rational_one;
    while (true)
    {
        new_generator->terms[a->generator->term_count] = number_rational_multiply(local_stack,
            output_stack, new_term, k);
        new_generator->minimal_polynomial = sum_minimal_polynomial(local_stack, output_stack,
            new_generator, a->generator->minimal_polynomial,
            new_generator->terms[a->generator->term_count]->minimal_polynomial);
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
        if (a->value.numerator->value_count == 0)
        {
            return number_copy(output_stack, b);
        }
        if (b->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*out = number_rational_initialize(output_stack,
                rational_add(local_stack, output_stack, &a->value, &b->value));
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
            out->term_count = 2;
            out->generator = number_copy(output_stack, a);
            out->terms = ARRAY_ALLOCATE(output_stack, 2, struct Number);
            out->terms[0] = out->generator;
            out->terms[1] = number_copy(output_stack, b);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, 2, struct RationalPolynomial);
            out->terms_in_terms_of_generator[0] = polynomial_allocate(output_stack, 2);
            out->terms_in_terms_of_generator[0]->coefficients[0] = &rational_zero;
            out->terms_in_terms_of_generator[0]->coefficients[1] = &rational_one;
            out->terms_in_terms_of_generator[1] = polynomial_allocate(output_stack, 1);
            out->terms_in_terms_of_generator[1]->coefficients[0] =
                rational_copy(output_stack, &b->value);
            if (!out->minimal_polynomial)
            {
                out->minimal_polynomial = sum_minimal_polynomial(output_stack, local_stack, out,
                    a->minimal_polynomial, b->minimal_polynomial);
            }
            return out;
        case '^':
        case '*':
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*a_term = ALLOCATE(local_stack, struct Number);
            a_term->operation = '+';
            a_term->generator = number_copy(local_stack, a);
            a_term->term_count = 1;
            a_term->terms = ALLOCATE(local_stack, struct Number*);
            a_term->terms[0] = a;
            a_term->terms_in_terms_of_generator =
                ALLOCATE(local_stack, struct RationalPolynomial*);
            struct { size_t coefficient_count; struct Rational*coefficients[2]; }x =
                { 2, { &rational_zero, &rational_one } };
            a_term->terms_in_terms_of_generator[0] = (struct RationalPolynomial*)&x;
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
            for (size_t i = b->term_count; i-- > 1;)
            {
                a = number_add(local_stack, output_stack, a, b->terms[i]);
            }
            a = number_add(output_stack, local_stack, a, b->terms[0]);
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

struct Number*number_rational_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Rational*b)
{
    if (b->numerator->value_count == 0)
    {
        return number_rational_initialize(output_stack, &rational_zero);
    }
    if (integer_equals(b->numerator, &one) && integer_equals(b->denominator, &one))
    {
        return number_copy(output_stack, a);
    }
    struct Number*out;
    switch (a->operation)
    {
    case 'r':
    {
        void*local_stack_savepoint = local_stack->cursor;
        out = number_rational_initialize(output_stack,
            rational_multiply(local_stack, output_stack, &a->value, b));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
        if (a->left->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            out = number_multiply(output_stack, local_stack,
                number_rational_multiply(local_stack, output_stack, a->left, b), a->right);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
    case '^':
    {
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '*';
        out->left = number_rational_initialize(output_stack, b);
        out->right = number_copy(output_stack, a);
        break;
    }
    default:
        out = ALLOCATE(output_stack, struct Number);
        out->operation = a->operation;
        if (a->operation == '+')
        {
            out->generator = number_copy(output_stack, a->generator);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, a->term_count, struct RationalPolynomial*);
            for (size_t i = 0; i < a->term_count; ++i)
            {
                out->terms_in_terms_of_generator[i] =
                    rational_polynomial_rational_multiply(output_stack, local_stack,
                        a->terms_in_terms_of_generator[i], b);
            }
        }
        out->term_count = a->term_count;
        out->terms = ARRAY_ALLOCATE(output_stack, a->term_count, struct Number*);
        for (size_t i = 0; i < a->term_count; ++i)
        {
            out->terms[i] = number_rational_multiply(output_stack, local_stack, a->terms[i], b);
        }
    }
    void*local_stack_savepoint = local_stack->cursor;
    out->minimal_polynomial =
        polynomial_allocate(output_stack, a->minimal_polynomial->coefficient_count);
    struct Rational*power = &rational_one;
    for (size_t i = 0; i < a->minimal_polynomial->coefficient_count; ++i)
    {
        out->minimal_polynomial->coefficients[i] = rational_divide(output_stack, local_stack,
            a->minimal_polynomial->coefficients[i], power);
        power = rational_unreduced_multiply(local_stack, output_stack, power, b, 0);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_consolidate_product(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Number*b)
{
    if (b->operation == 'r')
    {
        return number_rational_multiply(output_stack, local_stack, a, &b->value);
    }
    switch (a->operation)
    {
    case 'r':
    {
        return number_rational_multiply(output_stack, local_stack, b, &a->value);
    }
    case '^':
    {
        if (b->operation == '^')
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct ExtendedGCDInfo info;
            integer_extended_gcd(local_stack, output_stack, &info, a->right->value.denominator,
                b->right->value.denominator);
            struct Integer*product_index = integer_euclidean_quotient(local_stack, output_stack,
                integer_multiply(local_stack, output_stack, a->right->value.denominator,
                    b->right->value.denominator), info.gcd);
            struct Number*radicand_factor_a = number_exponentiate(local_stack, output_stack,
                a->left, &(struct Rational){info.b_over_gcd, &one});
            struct Number*radicand_factor_b = number_exponentiate(local_stack, output_stack,
                b->left, &(struct Rational){info.a_over_gcd, &one});
            struct RationalInterval radicand_factor_a_argument_estimate;
            number_rational_argument_estimate(local_stack, output_stack,
                &radicand_factor_a_argument_estimate, radicand_factor_a, &rational_one);
            struct RationalInterval radicand_factor_b_argument_estimate;
            number_rational_argument_estimate(local_stack, output_stack,
                &radicand_factor_b_argument_estimate, radicand_factor_b, &rational_one);
            if (!integer_equals(a->right->value.denominator, b->right->value.denominator))
            {
                struct Rational*radicand_factor_argument_interval_size =
                    rational_integer_divide(local_stack, output_stack, pi_estimate_min,
                        product_index);
                struct RationalInterval a_argument_estimate;
                number_rational_argument_estimate(local_stack, output_stack, &a_argument_estimate,
                    a, radicand_factor_argument_interval_size);
                struct RationalInterval b_argument_estimate;
                number_rational_argument_estimate(local_stack, output_stack, &b_argument_estimate,
                    b, radicand_factor_argument_interval_size);
                if (rational_compare(output_stack, local_stack,
                    rational_integer_divide(local_stack, output_stack,
                        radicand_factor_a_argument_estimate.max, product_index),
                    a_argument_estimate.min) < 0 ||
                    rational_compare(output_stack, local_stack,
                        rational_integer_divide(local_stack, output_stack,
                            radicand_factor_b_argument_estimate.max, product_index),
                        b_argument_estimate.min) < 0)
                {
                    local_stack->cursor = local_stack_savepoint;
                    return number_product_consolidation_failed;
                }
            }
            struct Number*product_radicand =
                number_multiply(local_stack, output_stack, radicand_factor_a, radicand_factor_b);
            struct Number*out = number_exponentiate(local_stack, output_stack, product_radicand,
                &(struct Rational){&one, product_index});
            struct RationalInterval product_radicand_rational_argument_estimate;
            number_rational_argument_estimate(local_stack, output_stack,
                &product_radicand_rational_argument_estimate, product_radicand, pi_estimate_min);
            if (rational_compare(output_stack, local_stack,
                product_radicand_rational_argument_estimate.max,
                rational_add(local_stack, output_stack, radicand_factor_a_argument_estimate.min,
                    radicand_factor_b_argument_estimate.min)) < 0)
            {
                struct Number*root_of_unity =
                    get_roots_of_unity(output_stack, local_stack, product_index)[1];
                out = number_multiply(output_stack, local_stack, out, root_of_unity);
            }
            else
            {
                out = number_copy(output_stack, out);
            }
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        break;
    }
    case '*':
        switch (b->operation)
        {
        case '^':
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*factor =
                number_consolidate_product(local_stack, output_stack, a->left, b);
            if (factor)
            {
                struct Number*out = number_multiply(output_stack, local_stack, a->right, factor);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
            factor = number_consolidate_product(local_stack, output_stack, a->right, b);
            if (factor)
            {
                struct Number*out = number_multiply(output_stack, local_stack, a->left, factor);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
            local_stack->cursor = local_stack_savepoint;
            return number_product_consolidation_failed;
        }
        case '*':
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*out = number_multiply(output_stack, local_stack, a->left,
                number_multiply(local_stack, output_stack, a->right, b));
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        }
        break;
    case '+':
    case 'g':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = b;
        for (size_t i = a->term_count; i-- > 1; ++i)
        {
            out = number_add(local_stack, output_stack, out,
                number_multiply(local_stack, output_stack, a->terms[i], b));
        }
        out = number_add(output_stack, local_stack, out,
            number_multiply(local_stack, output_stack, a->terms[0], b));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    }
    return number_multiply(output_stack, local_stack, b, a);
}

struct Number*number_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Number*b)
{
    struct Number*out = number_consolidate_product(output_stack, local_stack, a, b);
    if (out == number_product_consolidation_failed)
    {
        void*local_stack_savepoint = local_stack->cursor;
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '*';
        out->left = number_copy(output_stack, a);
        out->right = number_copy(output_stack, b);
        struct NestedPolynomial*t =
            polynomial_allocate(local_stack, a->minimal_polynomial->coefficient_count);
        for (size_t i = 0; i < a->minimal_polynomial->coefficient_count; ++i)
        {
            t->coefficients[i] =
                polynomial_allocate(local_stack, a->minimal_polynomial->coefficient_count - i);
            size_t degree = t->coefficients[i]->coefficient_count - 1;
            for (size_t j = 0; j < degree; ++j)
            {
                t->coefficients[i]->coefficients[j] = &rational_zero;
            }
            t->coefficients[i]->coefficients[degree] =
                a->minimal_polynomial->coefficients[degree - i];
        }
        out->minimal_polynomial =
            number_minimal_polynomial_from_annulling_polynomial(output_stack, local_stack,
                nested_polynomial_resultant(local_stack, output_stack, t,
                    rational_polynomial_to_nested_polynomial(local_stack, b->minimal_polynomial)),
                out);
        local_stack->cursor = local_stack_savepoint;
    }
    return out;
}

struct Number*number_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out;
    switch (a->operation)
    {
    case 'r':
    {
        if (!a->value.numerator->value_count)
        {
            puts("Tried to divide by 0.");
            return number_divide_by_zero_error;
        }
        out = number_rational_initialize(output_stack,
            rational_reciprocal(local_stack, output_stack, &a->value));
        break;
    }
    case '^':
    {
        out = number_divide(output_stack, local_stack,
            number_exponentiate(local_stack, output_stack, a->left,
                &(struct Rational){ integer_add(local_stack, a->right->value.denominator,
                    &INT(1, -)), a->right->value.denominator }), a->left);
        break;
    }
    case '*':
    {
        out = number_multiply(output_stack, local_stack,
            number_reciprocal(local_stack, output_stack, a->left),
            number_reciprocal(local_stack, output_stack, a->right));
        break;
    }
    case '+':
    {
        struct Number**conjugates = number_conjugates(local_stack, output_stack, a);
        out = &number_one;
        for (size_t i = 0; i < a->minimal_polynomial->coefficient_count - 1; ++i)
        {
            out = number_multiply(local_stack, output_stack, out, conjugates[i]);
        }
        struct Rational*coefficient =
            rational_reciprocal(local_stack, output_stack, a->minimal_polynomial->coefficients[0]);
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
    if (reciprocal == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    return number_multiply(output_stack, local_stack, dividend, reciprocal);
}

struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        return rational_copy(output_stack, &a->value);
    case '^':
        return &rational_one;
    case '*':
        return number_rational_factor(output_stack, local_stack, a->left);
    case '+':
    case 'g':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational*out = &rational_one;
        for (size_t i = 0; i < a->term_count; ++i)
        {
            struct Rational*term_factor =
                number_rational_factor(local_stack, output_stack, a->terms[i]);
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

struct Number*number_exponentiate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*base, struct Rational*exponent)
{
    if (!exponent->numerator->value_count)
    {
        return &number_one;
    }
    if (exponent->numerator->sign < 0)
    {
        void*local_stack_savepoint = local_stack->cursor;
        base = number_reciprocal(local_stack, output_stack, base);
        if (base == number_divide_by_zero_error)
        {
            local_stack->cursor = local_stack_savepoint;
            return number_divide_by_zero_error;
        }
        struct Number*out = number_exponentiate(output_stack, local_stack, base,
            rational_negative(local_stack, exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    switch (base->operation)
    {
    case 'r':
    {
        if (base->value.numerator->value_count == 0)
        {
            return number_copy(output_stack, base);
        }
        void*local_stack_savepoint = local_stack->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*out = number_rational_initialize(output_stack,
                rational_exponentiate(local_stack, output_stack, &base->value,
                    exponent->numerator));
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (!integer_equals(base->value.denominator, &one))
        {
            struct Number*out = number_rational_multiply(output_stack, local_stack,
                number_exponentiate(local_stack, output_stack,
                    number_rational_initialize(local_stack,
                        &(struct Rational){integer_multiply(local_stack, output_stack,
                            base->value.numerator,
                            integer_exponentiate(local_stack, output_stack, base->value.denominator,
                                integer_add(local_stack, exponent->denominator, &INT(1, -)))),
                        &one}), exponent), &(struct Rational){&one, base->value.denominator});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (!integer_equals(exponent->denominator, &one))
        {
            struct Integer*radicand = integer_magnitude(local_stack, base->value.numerator);
            struct Factor*factors;
            size_t factor_count = integer_factor(local_stack, output_stack, &factors, radicand);
            struct Integer*coefficient = &one;
            struct Integer*multiplicity_gcd = exponent->denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct IntegerDivision multiplicity_reduction;
                integer_euclidean_divide(local_stack, output_stack, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(local_stack, output_stack, coefficient,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = integer_gcd(local_stack, output_stack, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            struct Integer*new_degree;
            if (base->value.numerator->sign > 0)
            {
                new_degree = integer_euclidean_quotient(local_stack, output_stack,
                    exponent->denominator, multiplicity_gcd);
                if (integer_equals(new_degree, &one))
                {
                    struct Rational*rational_radicand = ALLOCATE(output_stack, struct Rational);
                    rational_radicand->numerator = integer_copy(output_stack, radicand);
                    rational_radicand->denominator = &one;
                    struct Number*out = number_rational_initialize(output_stack, rational_radicand);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
            else
            {
                new_degree = exponent->denominator;
            }
            radicand = integer_initialize(local_stack, 1, base->value.numerator->sign);
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
                    number_rational_initialize(local_stack, &(struct Rational){&one, new_degree})),
                &(struct Rational){coefficient, &one});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        return base;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_exponentiate(output_stack, local_stack, base->left,
            rational_multiply(local_stack, output_stack, exponent, &base->right->value));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_multiply(output_stack, local_stack,
            number_exponentiate(local_stack, output_stack, base->left, exponent),
            number_exponentiate(local_stack, output_stack, base->right, exponent));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*out = number_exponentiate(output_stack, local_stack,
                generic_exponentiate(&(struct RingOperations){number_copy, 0, 0, &number_one, 0,
                    0, number_generic_multiply},
                    local_stack, output_stack, base, exponent->numerator, 0),
                &(struct Rational){&one, exponent->denominator});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Rational*base_rational_factor =
            number_rational_factor(local_stack, output_stack, base);
        struct Factor*factors;
        size_t factor_count =
            integer_factor(local_stack, output_stack, &factors, base_rational_factor->numerator);
        struct Rational base_cancelling_rational_factor =
            { integer_exponentiate(local_stack, output_stack, base_rational_factor->denominator,
                exponent->denominator), &one };
        struct Rational product_rational_factor = { &one, base_rational_factor->denominator };
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_quotient(local_stack,
                output_stack, factors[factor_index].multiplicity, exponent->denominator);
            base_cancelling_rational_factor.denominator =
                integer_multiply(local_stack, output_stack,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        integer_multiply(local_stack, output_stack, reduced_multiplicity,
                            exponent->denominator)), base_cancelling_rational_factor.denominator);
            product_rational_factor.numerator =
                integer_multiply(local_stack, output_stack, product_rational_factor.numerator,
                    integer_exponentiate(local_stack, output_stack, factors[factor_index].value,
                        reduced_multiplicity));
        }
        struct Number*out = number_rational_multiply(output_stack, local_stack,
            number_surd_initialize(local_stack, output_stack,
                number_rational_multiply(local_stack, output_stack, base,
                    &base_cancelling_rational_factor), 
                number_rational_initialize(local_stack, exponent)), &product_rational_factor);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number**number_sum_or_product_conjugates(struct Number*(operation)(struct Stack*,
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
            struct Number*candidate =
                operation(local_stack, output_stack, a_conjugates[i], b_conjugates[j]);
            if (rational_polynomial_equals(minimal_polynomial, candidate->minimal_polynomial))
            {
                out[conjugate_count] = number_copy(output_stack, candidate);
                ++conjugate_count;
                if (conjugate_count == minimal_polynomial->coefficient_count - 1)
                {
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
        }
    }
    crash("Not enough conjugates found.");
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
        struct Number**radicand_conjugates = number_conjugates(local_stack, output_stack, a->left);
        for (size_t i = 0; i < a->left->minimal_polynomial->coefficient_count - 1; ++i)
        {
            radicand_conjugates[i] = number_exponentiate(local_stack, output_stack,
                radicand_conjugates[i], &a->right->value);
        }
        size_t roots_of_unity_count = integer_to_size_t(a->right->value.denominator);
        struct Number**roots_of_unity =
            get_roots_of_unity(output_stack, local_stack, a->right->value.denominator);
        struct Number**out = number_sum_or_product_conjugates(number_multiply, output_stack,
            local_stack, radicand_conjugates, a->left->minimal_polynomial->coefficient_count - 1,
            roots_of_unity, roots_of_unity_count, a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number**left_conjugates = number_conjugates(local_stack, output_stack, a->left);
        struct Number**right_conjugates = number_conjugates(local_stack, output_stack, a->right);
        struct Number**out = number_sum_or_product_conjugates(number_multiply, output_stack,
            local_stack, left_conjugates, a->left->minimal_polynomial->coefficient_count - 1,
            right_conjugates, a->right->minimal_polynomial->coefficient_count - 1,
            a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number**out = number_conjugates(local_stack, output_stack, a->terms[0]);
        struct RationalPolynomial*left_minimal_polynomial = a->terms[0]->minimal_polynomial;
        size_t original_term_count = a->term_count;
        a->term_count = 2;
        while (a->term_count < original_term_count)
        {
            struct Number*right_term = a->terms[a->term_count - 1];
            struct RationalPolynomial*new_left_minimal_polynomial =
                sum_minimal_polynomial(local_stack, output_stack, a, left_minimal_polynomial,
                    right_term->minimal_polynomial);
            out = number_sum_or_product_conjugates(number_add, local_stack, output_stack, out,
                left_minimal_polynomial->coefficient_count - 1,
                number_conjugates(local_stack, output_stack, right_term),
                right_term->minimal_polynomial->coefficient_count - 1, new_left_minimal_polynomial);
            left_minimal_polynomial = new_left_minimal_polynomial;
            ++a->term_count;
        }
        out = number_sum_or_product_conjugates(number_add, output_stack, local_stack, out,
            left_minimal_polynomial->coefficient_count - 1,
            number_conjugates(local_stack, output_stack, a->terms[a->term_count - 1]),
            a->terms[a->term_count - 1]->minimal_polynomial->coefficient_count - 1,
            a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_evaluate(struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    if (a->operation == 'r')
    {
        a->minimal_polynomial = rational_minimal_polynomial(output_stack, &a->value);
        return a;
    }
    void*local_stack_savepoint = local_stack->cursor;
    a->left = number_evaluate(local_stack, output_stack, a->left);
    if (a->left == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    a->right = number_evaluate(local_stack, output_stack, a->right);
    if (a->right == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    struct Number*out;
    switch (a->operation)
    {
    case '+':
        out = number_add(output_stack, local_stack, a->left, a->right);
        break;
    case '-':
    {
        out = number_add(output_stack, local_stack, a->left,
            number_rational_multiply(local_stack, output_stack, a->right,
                &(struct Rational){ integer_initialize(local_stack, 1, -1), &one }));
        break;
    }
    case '*':
        out = number_multiply(output_stack, local_stack, a->left, a->right);
        break;
    case '/':
        out = number_divide(output_stack, local_stack, a->left, a->right);
        break;
    case '^':
        if (a->right->operation != 'r')
        {
            puts("The input expression contains an exponentiation whose exponent is not both real "
                "and rational; this program doesn't handle transcendental numbers.");
            return 0;
        }
        out = number_exponentiate(output_stack, local_stack, a->left, &a->right->value);
        break;
    default:
        crash("Number operation not recognized.");
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

size_t number_string(struct Stack*output_stack, struct Stack*local_stack, struct Number*number)
{
    switch (number->operation)
    {
    case 'r':
    {
        size_t char_count = integer_string(output_stack, local_stack, number->value.numerator);
        if (!integer_equals(number->value.denominator, &one))
        {
            *(char*)ALLOCATE(output_stack, char) = '/';
            char_count += 1 + integer_string(output_stack, local_stack, number->value.denominator);
        }
        return char_count;
    }
    case '^':
    {
        size_t char_count;
        if (number->left->operation != 'r' || number->left->value.numerator->sign < 0)
        {
            *(char*)ALLOCATE(output_stack, char) = '(';
            char_count = 2 + number_string(output_stack, local_stack, number->left);
            *(char*)ALLOCATE(output_stack, char) = ')';
        }
        else
        {
            char_count = number_string(output_stack, local_stack, number->left);
        }
        char*exponent = stack_slot_allocate(output_stack, 2 * sizeof(char), _Alignof(char));
        exponent[0] = '^';
        exponent[1] = '(';
        char_count += 3 + number_string(output_stack, local_stack, number->right);
        *(char*)ALLOCATE(output_stack, char) = ')';
        return char_count;
    }
    case '*':
    {
        size_t char_count;
        if (number->left->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            if (integer_equals(number->left->value.numerator,
                integer_initialize(local_stack, 1, -1)))
            {
                *(char*)ALLOCATE(output_stack, char) = '-';
                char_count = 1 + number_string(output_stack, local_stack, number->right);
            }
            else if (!integer_equals(number->left->value.numerator, &one))
            {
                char_count =
                    integer_string(output_stack, local_stack, number->left->value.numerator);
                char*right_string = output_stack->cursor;
                size_t right_char_count = number_string(output_stack, local_stack, number->right);
                if (right_string[0] != '(')
                {
                    ALLOCATE(output_stack, char);
                    memcpy(right_string + 1, right_string, right_char_count);
                    *right_string = '*';
                    char_count += 1 + right_char_count;
                }
            }
            else
            {
                char_count = number_string(output_stack, local_stack, number->right);
            }
            if (!integer_equals(number->left->value.denominator, &one))
            {
                *(char*)ALLOCATE(output_stack, char) = '/';
                char_count +=
                    1 + integer_string(output_stack, local_stack, number->left->value.denominator);
            }
            local_stack->cursor = local_stack_savepoint;
        }
        else
        {
            char_count = number_string(output_stack, local_stack, number->left) + 1;
            *(char*)ALLOCATE(output_stack, char) = '*';
            char_count += number_string(output_stack, local_stack, number->right);
        }
        return char_count;
    }
    case '+':
    {
        size_t char_count = number_string(output_stack, local_stack, number->left) + 1;
        *(char*)ALLOCATE(output_stack, char) = '+';
        return char_count + number_string(output_stack, local_stack, number->right);
    }
    default:
        crash("Number operation not recognized.");
    }
}