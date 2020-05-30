#include "declarations.h"

struct TermSplit
{
    struct Rational rational_part;
    struct Number*irrational_part;
};

struct RationalPolynomial*sum_get_minimal_polynomial(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a,
    struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*power = nested_polynomial_one;
    POLY(d_term0, struct RationalPolynomial, struct Rational, 2, rational_zero, rational_one);
    POLY(d_term1, struct RationalPolynomial, struct Rational, 1, { INT(1, -1), &one });
    POLY(d, struct NestedPolynomial, struct RationalPolynomial*, 2, &d_term0.p, &d_term1.p);
    struct NestedPolynomial*t = polynomial_zero;
    for (size_t i = 0; i < left_minimal_polynomial->coefficient_count; ++i)
    {
        POLY(factor, struct RationalPolynomial, struct Rational, 1,
            left_minimal_polynomial->coefficients[i]);
        t = nested_polynomial_add(local_stack, output_stack, t,
            nested_polynomial_rational_polynomial_multiply(local_stack, output_stack, power,
                &factor.p));
        power = nested_polynomial_multiply(local_stack, output_stack, power, &d.p);
    }
    struct RationalPolynomial*out =
        number_annulling_polynomial_to_minimal_polynomial(output_stack, local_stack,
            nested_polynomial_get_resultant(local_stack, output_stack, t,
                rational_polynomial_to_nested_polynomial(local_stack, right_minimal_polynomial)),
            a);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void term_split(struct Stack*output_stack, struct TermSplit*out, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        out->rational_part = rational_copy(output_stack, &a->value);
        out->irrational_part = &number_one;
        return;
    case '^':
        out->rational_part = rational_one;
        out->irrational_part = number_copy(output_stack, a);
        return;
    case '*':
        if (a->elements[0]->operation == 'r')
        {
            out->rational_part = rational_copy(output_stack, &a->elements[0]->value);
            if (a->element_count == 2)
            {
                out->irrational_part = number_copy(output_stack, a->elements[1]);
            }
            else
            {
                out->irrational_part = number_copy(output_stack, a);
                out->irrational_part->elements += 1;
                out->irrational_part->element_count -= 1;
            }
        }
        else
        {
            out->irrational_part = number_copy(output_stack, a);
        }
    }
}

struct Number*term_consolidate(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Number*b)
{
    if (a->operation == 'r' && !a->value.numerator->value_count)
    {
        return number_copy(output_stack, b);
    }
    if (b->operation == 'r' && !b->value.numerator->value_count)
    {
        return number_copy(output_stack, a);
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct TermSplit a_split;
    term_split(local_stack, &a_split, a);
    struct TermSplit b_split;
    term_split(local_stack, &b_split, b);
    if (!number_formal_compare(output_stack, local_stack, a_split.irrational_part,
        b_split.irrational_part))
    {
        struct Rational sum =
            rational_add(local_stack, output_stack, &a_split.rational_part, &b_split.rational_part);
        struct Number*out =
            number_rational_multiply(output_stack, local_stack, a_split.irrational_part, &sum);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    local_stack->cursor = local_stack_savepoint;
    return 0;
}

struct Number*number_add(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a, struct Number*b)
{
    if (a->operation == '+')
    {
        void*local_stack_savepoint = local_stack->cursor;
        if (b->operation == '+')
        {
            struct Number*out = a;
            for (size_t i = b->element_count; i-- > 1;)
            {
                out = number_add(local_stack, output_stack, out, b->elements[i]);
            }
            out = number_add(output_stack, local_stack, out, b->elements[0]);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        struct Number**old_terms = ARRAY_ALLOCATE(local_stack, a->element_count, struct Number*);
        memcpy(old_terms, a->elements, a->element_count * sizeof(struct Number*));
        size_t term_count = a->element_count + 1;
        struct Number*consolidation = b;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            struct Number*consolidation_attempt =
                term_consolidate(local_stack, output_stack, a->elements[i], b);
            if (consolidation_attempt)
            {
                if (consolidation_attempt->operation == 'r' &&
                    !consolidation_attempt->value.numerator->value_count)
                {
                    local_stack->cursor = local_stack_savepoint;
                    if (a->element_count == 2)
                    {
                        return number_copy(output_stack, a->elements[1 - i]);
                    }
                    struct Number*out = ALLOCATE(output_stack, struct Number);
                    out->operation = '+';
                    out->element_count = a->element_count - 1;
                    out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number);
                    for (size_t j = 0; j < i; ++j)
                    {
                        out->elements[j] = number_copy(output_stack, a->elements[j]);
                    }
                    for (size_t j = i + 1; j < a->element_count; ++j)
                    {
                        out->elements[j - 1] = number_copy(output_stack, a->elements[j]);
                    }
                    out->minimal_polynomial = 0;
                    out->generator = 0;
                    return out;
                }
                old_terms[i] = 0;
                term_count = a->element_count;
                consolidation = consolidation_attempt;
                break;
            }
        }
        struct Number*out = ALLOCATE(output_stack, struct Number);
        out->operation = '+';
        out->elements = ARRAY_ALLOCATE(output_stack, term_count, struct Number);
        out->element_count = 0;
        out->minimal_polynomial = 0;
        out->generator = 0;
        for (size_t i = 0; i < a->element_count; ++i)
        {
            if (old_terms[i])
            {
                if (number_formal_compare(output_stack, local_stack, consolidation,
                    old_terms[i]) > 0)
                {
                    out->elements[out->element_count] = number_copy(output_stack, consolidation);
                    ++out->element_count;
                    for (size_t j = i; j < a->element_count; ++j)
                    {
                        if (old_terms[j])
                        {
                            out->elements[out->element_count] =
                                number_copy(output_stack, old_terms[j]);
                            ++out->element_count;
                        }
                    }
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
                else
                {
                    out->elements[out->element_count] = number_copy(output_stack, old_terms[i]);
                    ++out->element_count;
                }
            }
        }
        out->elements[out->element_count] = number_copy(output_stack, consolidation);
        ++out->element_count;
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    if (b->operation != '+')
    {
        struct Number*out = term_consolidate(output_stack, local_stack, a, b);
        if (!out)
        {
            out = ALLOCATE(output_stack, struct Number);
            out->operation = '+';
            out->element_count = 2;
            out->elements = ARRAY_ALLOCATE(output_stack, 2, struct Number);
            out->elements[0] = number_copy(output_stack, a);
            out->elements[1] = number_copy(output_stack, b);
            if (number_formal_compare(output_stack, local_stack, a, b) < 0)
            {
                SWAP(out->elements[0], out->elements[1], struct Number*);
            }
            out->minimal_polynomial = 0;
            out->generator = 0;
        }
        return out;
    }
    return number_add(output_stack, local_stack, b, a);
}

struct RationalPolynomial**sum_convert_terms_to_new_generator(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, size_t term_count,
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
        if (terms_in_terms_of_old_generator[i]->coefficient_count > max_term_coefficient_count)
        {
            max_term_coefficient_count = terms_in_terms_of_old_generator[i]->coefficient_count;
        }
    }
    struct RationalPolynomial**old_generator_powers =
        ARRAY_ALLOCATE(local_stack, max_term_coefficient_count, struct RationalPolynomial*);
    old_generator_powers[0] = rational_polynomial_one;
    for (size_t i = 1; i < max_term_coefficient_count; ++i)
    {
        old_generator_powers[i] = number_field_element_multiply(local_stack, output_stack,
            old_generator_powers[i - 1], old_generator_in_terms_of_new,
            new_generator_minimal_polynomial);
    }
    for (size_t i = 0; i < term_count; ++i)
    {
        out[i] = polynomial_zero;
        for (size_t j = 0; j < terms_in_terms_of_old_generator[i]->coefficient_count; ++j)
        {
            out[i] = rational_polynomial_add(local_stack, output_stack,
                rational_polynomial_rational_multiply(local_stack, output_stack,
                    old_generator_powers[j], terms_in_terms_of_old_generator[i]->coefficients + j),
                out[i]);
        }
        out[i] = rational_polynomial_copy(output_stack, out[i]);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_get_a_in_terms_of_b(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a, struct Number*b)
{
    if ((b->minimal_polynomial->coefficient_count - 1) %
        (a->minimal_polynomial->coefficient_count - 1) != 0)
    {
        return 0;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*nested_a_minimal_polynomial =
        rational_polynomial_to_nested_polynomial(local_stack,
            rational_polynomial_euclidean_divide(local_stack, output_stack, a->minimal_polynomial,
                b->minimal_polynomial).remainder);
    switch (nested_a_minimal_polynomial->coefficient_count)
    {
    case 0:
    {
        struct RationalPolynomial*out = POLYNOMIAL_ALLOCATE(output_stack, 2, struct Rational);
        out->coefficients[0] = rational_zero;
        out->coefficients[1] = rational_one;
        return out;
    }
    case 1:
        local_stack->cursor = local_stack_savepoint;
        return 0;
    }
    struct NestedPolynomial**a_minimal_polynomial_factors = ARRAY_ALLOCATE(local_stack,
        nested_a_minimal_polynomial->coefficient_count - 1, struct NestedPolynomial*);
    size_t candidate_factor_count = number_field_polynomial_factor(local_stack, output_stack,
        a_minimal_polynomial_factors, nested_a_minimal_polynomial, b->minimal_polynomial);
    struct RationalPolynomial**candidate_factors =
        (struct RationalPolynomial**)a_minimal_polynomial_factors;
    for (size_t i = 0; i < candidate_factor_count;)
    {
        if (a_minimal_polynomial_factors[i]->coefficient_count == 2)
        {
            candidate_factors[i] = rational_polynomial_negate(local_stack,
                a_minimal_polynomial_factors[i]->coefficients[0]);
            ++i;
        }
        else
        {
            --candidate_factor_count;
            a_minimal_polynomial_factors[i] = a_minimal_polynomial_factors[candidate_factor_count];
        }
    }
    struct Region a_estimate =
        number_get_region_estimate(local_stack, output_stack, a, &rational_one);
    while (rational_polynomial_count_roots_in_region(local_stack, output_stack,
        a->minimal_polynomial, &a_estimate) > 1)
    {
        number_halve_region_estimate_dimensions(local_stack, output_stack, &a_estimate, a);
    }
    struct Region b_estimate =
        number_get_region_estimate(local_stack, output_stack, b, &rational_one);
    for (size_t i = 0; i < candidate_factor_count; ++i)
    {
        struct Rational interval_size_for_evaluation =
        { &one, size_t_to_integer(local_stack, candidate_factors[i]->coefficient_count) };
        struct Region factor_at_b;
        while (true)
        {
            factor_at_b = rational_polynomial_evaluate_at_region_estimate(local_stack, output_stack,
                candidate_factors[i], &b_estimate, &interval_size_for_evaluation);
            if (regions_are_disjoint(output_stack, local_stack, &a_estimate, &factor_at_b))
            {
                goto candidate_rejected;
            }
            if (rational_polynomial_count_roots_in_region(local_stack, output_stack,
                a->minimal_polynomial, &factor_at_b) == 1)
            {
                if (region_a_contains_b(output_stack, local_stack, &a_estimate, &factor_at_b))
                {
                    struct RationalPolynomial*out =
                        rational_polynomial_copy(output_stack, candidate_factors[i]);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
                break;
            }
            number_halve_region_estimate_dimensions(local_stack, output_stack, &b_estimate, b);
            interval_size_for_evaluation.denominator =
                integer_double(local_stack, interval_size_for_evaluation.denominator);
        }
        while (true)
        {
            if (region_a_contains_b(output_stack, local_stack, &factor_at_b, &a_estimate))
            {
                struct RationalPolynomial*out =
                    rational_polynomial_copy(output_stack, candidate_factors[i]);
                local_stack->cursor = local_stack_savepoint;
                return out;
            }
            number_halve_region_estimate_dimensions(local_stack, output_stack, &a_estimate, a);
            if (regions_are_disjoint(output_stack, local_stack, &a_estimate, &factor_at_b))
            {
                goto candidate_rejected;
            }
        }
    candidate_rejected:;
    }
    local_stack->cursor = local_stack_savepoint;
    return 0;
}

struct Number*sum_append_term(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
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
        if (number_formal_compare(output_stack, local_stack, out->elements[i - 1],
            out->elements[i]) >= 0)
        {
            return out;
        }
        SWAP(out->elements[i], out->elements[i - 1], struct Number*);
        SWAP(out->terms_in_terms_of_generator[i], out->terms_in_terms_of_generator[i - 1],
            struct RationalPolynomial*);
    }
    return out;
}

struct Number*sum_incorporate_term_in_terms_of_new_generator(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number**a_terms, size_t a_term_count,
    struct RationalPolynomial*a_minimal_polynomial, struct Number*new_term,
    struct Number*new_generator, struct RationalPolynomial**a_terms_in_terms_of_new_generator,
    struct RationalPolynomial*new_term_in_terms_of_new_generator)
{
    void*local_stack_savepoint = output_stack->cursor;
    size_t max_term_coefficient_count = 0;
    for (size_t i = 0; i < a_term_count; ++i)
    {
        if (a_terms_in_terms_of_new_generator[i]->coefficient_count > max_term_coefficient_count)
        {
            max_term_coefficient_count = a_terms_in_terms_of_new_generator[i]->coefficient_count;
        }
    }
    struct Matrix matrix =
    { ARRAY_ALLOCATE(local_stack, max_term_coefficient_count, struct Rational*),
        a_term_count, max_term_coefficient_count };
    for (size_t i = 0; i < matrix.height; ++i)
    {
        matrix.rows[i] = ARRAY_ALLOCATE(local_stack, matrix.width, struct Rational);
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
            matrix.rows[j][i] = rational_zero;
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
    struct Rational*augmentation = ARRAY_ALLOCATE(local_stack, matrix.height, struct Rational);
    memcpy(augmentation, new_term_in_terms_of_new_generator->coefficients,
        new_term_in_terms_of_new_generator->coefficient_count * sizeof(struct Rational));
    for (size_t i = new_term_in_terms_of_new_generator->coefficient_count; i < matrix.height; ++i)
    {
        augmentation[i] = rational_zero;
    }
    matrix_make_row_echelon_form(local_stack, output_stack, &matrix, augmentation);
    for (size_t i = matrix.height; i-- > matrix.width;)
    {
        if (augmentation[i].numerator->value_count != 0)
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
                terms_remaining * sizeof(struct Rational));
            memmove(matrix.rows + i, matrix.rows + i + 1,
                terms_remaining * sizeof(struct Rational*));
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
        struct Rational scale =
            rational_add(local_stack, output_stack, augmentation + i, &rational_one);
        if (scale.numerator->value_count)
        {
            out->elements[out->element_count] =
                number_rational_multiply(output_stack, local_stack, a_terms[i], &scale);
            out->terms_in_terms_of_generator[out->element_count] =
                rational_polynomial_rational_multiply(output_stack, local_stack,
                    a_terms_in_terms_of_new_generator[i], &scale);
            ++out->element_count;
        }
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

//Can handle the degenerate case of a_term_count == 1, and returns a sum even in degenerate cases of
//fewer than two terms. Leaves it to the caller to fill out the minimal_polynomial field since it
//will need to be calculated only in nondegenerate cases.
struct Number*sum_incorporate_term(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number**a_terms, size_t a_term_count,
    struct RationalPolynomial**a_terms_in_terms_of_generator, struct Number*a_generator,
    struct RationalPolynomial*a_minimal_polynomial, struct Number*new_term)
{
    if (new_term->operation == 'r')
    {
        POLY(new_term_polynomial, struct RationalPolynomial, struct Rational, 1, new_term->value);
        return sum_append_term(output_stack, local_stack, a_terms, a_term_count,
            a_terms_in_terms_of_generator, new_term, a_generator, &new_term_polynomial.p);
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*new_term_in_terms_of_a_generator =
        number_get_a_in_terms_of_b(local_stack, output_stack, new_term, a_generator);
    if (new_term_in_terms_of_a_generator)
    {
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a_terms, a_term_count, a_minimal_polynomial, new_term, a_generator,
            a_terms_in_terms_of_generator, new_term_in_terms_of_a_generator);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    POLY(x, struct RationalPolynomial, struct Rational, 2, rational_zero, rational_one);
    struct RationalPolynomial*old_generator_in_terms_of_new_term =
        number_get_a_in_terms_of_b(local_stack, output_stack, a_generator, new_term);
    if (old_generator_in_terms_of_new_term)
    {
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
            local_stack, a_terms, a_term_count, a_minimal_polynomial, new_term, new_term,
            sum_convert_terms_to_new_generator(local_stack, output_stack, a_term_count,
                a_terms_in_terms_of_generator, old_generator_in_terms_of_new_term,
                new_term->minimal_polynomial), &x.p);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    struct Number*new_generator = ALLOCATE(local_stack, struct Number);
    new_generator->operation = '+';
    new_generator->generator = 0;
    if (a_generator->operation == '+')
    {
        new_generator->element_count = a_generator->element_count + 1;
        new_generator->elements =
            ARRAY_ALLOCATE(local_stack, new_generator->element_count, struct Number*);
        memcpy(new_generator->elements, a_generator->elements,
            a_generator->element_count * sizeof(struct Number*));
    }
    else
    {
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
        new_generator->minimal_polynomial = sum_get_minimal_polynomial(local_stack, output_stack,
            new_generator, a_generator->minimal_polynomial,
            new_generator->elements[new_generator->element_count - 1]->minimal_polynomial);
        struct RationalPolynomial*term_in_terms_of_new_generator =
            number_get_a_in_terms_of_b(local_stack, output_stack, new_term, new_generator);
        if (term_in_terms_of_new_generator)
        {
            struct Number*out = sum_incorporate_term_in_terms_of_new_generator(output_stack,
                local_stack, a_terms, a_term_count, a_minimal_polynomial, new_term, new_generator,
                sum_convert_terms_to_new_generator(local_stack, output_stack, a_term_count,
                    a_terms_in_terms_of_generator,
                    rational_polynomial_subtract(local_stack, output_stack, &x.p,
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

struct Number*number_incorporate_term(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number**a_terms, size_t a_term_count,
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
    out->minimal_polynomial = sum_get_minimal_polynomial(output_stack, local_stack, out,
        a_minimal_polynomial, new_term->minimal_polynomial);
    return out;
}

struct Number*number_eliminate_linear_dependencies(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a)
{
    if (a->minimal_polynomial)
    {
        return number_copy(output_stack, a);
    }
    void*local_stack_savepoint = local_stack->cursor;
    POLY(x_value, struct RationalPolynomial, struct Rational, 2, rational_zero, rational_one);
    struct RationalPolynomial*x = &x_value.p;
    struct Number*out = a->elements[0];
    for (size_t i = 1; i < a->element_count; ++i)
    {
        if (out->operation == '+')
        {
            out = number_incorporate_term(local_stack, output_stack, out->elements,
                out->element_count, out->terms_in_terms_of_generator, out->generator,
                out->minimal_polynomial, a->elements[i]);
        }
        else
        {
            out = number_incorporate_term(local_stack, output_stack, &out, 1, &x, out,
                out->minimal_polynomial, a->elements[i]);
        }
    }
    out = number_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}