#include "declarations.h"

struct Number*sum_incorporate_term_in_terms_of_new_generator(struct PoolSet*pool_set,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Number*new_term,
    struct Number*new_generator, struct Matrix*a_terms_in_terms_of_generator,
    struct RationalPolynomial*new_term_in_terms_of_generator)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Number*out = number_allocate(pool_set);
    out->operation = '+';
    out->generator = new_generator;
    increment_reference_count(a->first_term);
    out->first_term = a->first_term;
    struct Rational**augmentation =
        STACK_ARRAY_ALLOCATE(stack_a, a_terms_in_terms_of_generator->height, struct Rational*);
    memcpy(augmentation, new_term_in_terms_of_generator->coefficients,
        new_term_in_terms_of_generator->coefficient_count * sizeof(struct Rational*));
    for (size_t i = new_term_in_terms_of_generator->coefficient_count;
        i < a_terms_in_terms_of_generator->height; ++i)
    {
        augmentation[i] = &rational_zero;
    }
    matrix_row_echelon_form(rational_multiply, rational_subtract, stack_a, stack_b,
        a_terms_in_terms_of_generator, augmentation);
    for (size_t i = a_terms_in_terms_of_generator->height;
        i-- > a_terms_in_terms_of_generator->width;)
    {
        if (augmentation[i]->numerator->value_count != 0)
        {
            struct Term*term = out->first_term;
            while (term->next)
            {
                term = term->next;
            }
            term->next = pool_value_allocate(pool_set, sizeof(struct Term));
            term->next->value = new_term;
            term->next->in_terms_of_generator =
                rational_polynomial_copy_to_pool(pool_set, new_term_in_terms_of_generator);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        else
        {
            size_t terms_remaining = a_terms_in_terms_of_generator->height - i - 1;
            memcpy(augmentation + i, augmentation + i + 1, terms_remaining);
            memcpy(a_terms_in_terms_of_generator->rows + i,
                a_terms_in_terms_of_generator->rows + i + 1, terms_remaining);
            --a_terms_in_terms_of_generator->height;
        }
    }
    matrix_diagonalize(stack_a, stack_b, a_terms_in_terms_of_generator, augmentation);
    struct Term*previous_term = 0;
    struct Term*term = out->first_term;
    size_t out_term_count = 0;
    for (size_t i = 0; i < a_terms_in_terms_of_generator->width; ++i)
    {
        struct Rational*scale = rational_add(stack_a, stack_b, augmentation[i], &rational_one);
        if (scale->numerator->value_count)
        {
            term->value = number_rational_multiply(pool_set, stack_a, stack_b, term->value, scale);
            struct RationalPolynomial*in_terms_of_generator =
                rational_polynomial_rational_multiply(stack_a, stack_b, term->in_terms_of_generator,
                    scale);
            rational_polynomial_free(pool_set, term->in_terms_of_generator);
            term->in_terms_of_generator =
                rational_polynomial_copy_to_pool(pool_set, term->in_terms_of_generator);
            term = term->next;
            ++out_term_count;
        }
        else
        {
            if (previous_term)
            {
                previous_term->next = term->next;
                term_free(pool_set, term);
                term = previous_term->next;
            }
            else
            {
                previous_term = term->next;
                term_free(pool_set, term);
                term = previous_term;
            }
        }
    }
    switch (out_term_count)
    {
    case 0:
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value.numerator = &zero;
        out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
        number_free(pool_set, out);
        return out;
    }
    case 1:
    {
        struct Number*out_term = out->first_term->value;
        increment_reference_count(out_term);
        number_free(pool_set, out);
        return out_term;
    }
    default:
        return out;
    }
}

struct Term*sum_term_count_and_max_degree(struct Number*a, size_t*out_term_count,
    size_t*out_max_degree)
{
    *out_term_count = 0;
    *out_max_degree = 0;
    struct Term*term = a->first_term;
    while (true)
    {
        ++(*out_term_count);
        *out_max_degree = max(*out_max_degree, term->in_terms_of_generator->coefficient_count);
        if (!term->next)
        {
            return term;
        }
        term = term->next;
    }
}

void sun_convert_terms_to_new_generator(struct Stack*output_stack, struct Stack*local_stack,
    struct Matrix*out, struct Number*a, struct RationalPolynomial*old_generator_in_terms_of_new,
    struct RationalPolynomial*new_generator_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t term_max_degree = 0;
    sum_term_count_and_max_degree(a, &out->width, &term_max_degree);
    out->height = 0;
    struct RationalPolynomial**old_generator_powers =
        STACK_ARRAY_ALLOCATE(local_stack, out->height, struct RationalPolynomial*);
    old_generator_powers[0] = &rational_polynomial_one;
    for (size_t i = 1; i < term_max_degree; ++i)
    {
        old_generator_powers[i] = number_field_element_multiply(local_stack, output_stack,
            old_generator_powers[i - 1], old_generator_in_terms_of_new,
            new_generator_minimal_polynomial);
        out->height = max(out->height, old_generator_powers[i]->coefficient_count);
    }
    out->rows = STACK_ARRAY_ALLOCATE(output_stack, out->height, struct Rational**);
    for (size_t i = 0; i < out->height; ++i)
    {
        out->rows[i] = STACK_ARRAY_ALLOCATE(output_stack, out->width, struct Rational*);
    }
    struct Term*a_term = a->first_term;
    for (size_t i = 0; i < out->width; ++i)
    {
        struct RationalPolynomial*term_in_terms_of_new_generator =
            (struct RationalPolynomial*)&polynomial_zero;
        for (size_t j = 0; j < a_term->in_terms_of_generator->coefficient_count; ++j)
        {
            term_in_terms_of_new_generator = rational_polynomial_add(local_stack, output_stack,
                rational_polynomial_rational_multiply(local_stack, output_stack,
                    old_generator_powers[j], a_term->in_terms_of_generator->coefficients[j]),
                term_in_terms_of_new_generator);
        }
        for (size_t j = 0; j < term_in_terms_of_new_generator->coefficient_count; ++j)
        {
            out->rows[j][i] = rational_copy_to_stack(output_stack,
                term_in_terms_of_new_generator->coefficients[j]);
        }
        for (size_t j = term_in_terms_of_new_generator->coefficient_count; j < out->height; ++j)
        {
            out->rows[j][i] = &rational_zero;
        }
        a_term = a_term->next;
    }
}

struct Number*sum_incorporate_term(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Number*new_term)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct RationalPolynomial*new_term_in_terms_of_generator =
        number_a_in_terms_of_b(pool_set, stack_a, stack_b, new_term, a->generator);
    if (new_term_in_terms_of_generator)
    {
        struct Matrix a_terms_in_terms_of_generator;
        struct Term*a_last_term = sum_term_count_and_max_degree(a,
            &a_terms_in_terms_of_generator.width, &a_terms_in_terms_of_generator.height);
        if (new_term_in_terms_of_generator->coefficient_count >
            a_terms_in_terms_of_generator.height)
        {
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->generator = a->generator;
            increment_reference_count(a->generator);
            out->first_term = a->first_term;
            a->first_term = 0;
            number_free(pool_set, a);
            a_last_term->next = pool_value_allocate(pool_set, sizeof(struct Term));
            a_last_term->next->value = new_term;
            a_last_term->next->in_terms_of_generator =
                rational_polynomial_copy_to_pool(pool_set, new_term_in_terms_of_generator);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        a_terms_in_terms_of_generator.rows =
            STACK_ARRAY_ALLOCATE(stack_a, a_terms_in_terms_of_generator.height, struct Rational**);
        for (size_t i = 0; i < a_terms_in_terms_of_generator.height; ++i)
        {
            a_terms_in_terms_of_generator.rows[i] = STACK_ARRAY_ALLOCATE(stack_a,
                a_terms_in_terms_of_generator.width, struct Rational*);
        }
        struct Term*a_term = a->first_term;
        for (size_t i = 0; i < a_terms_in_terms_of_generator.width; ++i)
        {
            for (size_t j = 0; j < a_term->in_terms_of_generator->coefficient_count; ++j)
            {
                a_terms_in_terms_of_generator.rows[j][i] =
                    a_term->in_terms_of_generator->coefficients[j];
            }
            for (size_t j = a_term->in_terms_of_generator->coefficient_count;
                j < a_terms_in_terms_of_generator.height; ++j)
            {
                a_terms_in_terms_of_generator.rows[j][i] = &rational_zero;
            }
            a_term = a_term->next;
        }
        increment_reference_count(a->generator);
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(pool_set, stack_a,
            stack_b, a, new_term, a->generator, &a_terms_in_terms_of_generator,
            new_term_in_terms_of_generator);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    struct RationalPolynomial*x = stack_polynomial_allocate(stack_a, 2);
    x->coefficients[0] = &rational_zero;
    x->coefficients[0] = &rational_one;
    struct RationalPolynomial*old_generator_in_terms_of_new_term =
        number_a_in_terms_of_b(pool_set, stack_a, stack_b, a->generator, new_term);
    if (old_generator_in_terms_of_new_term)
    {
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, new_term);
        struct Matrix a_terms_in_terms_of_new_term;
        sun_convert_terms_to_new_generator(stack_a, stack_b, &a_terms_in_terms_of_new_term, a,
            old_generator_in_terms_of_new_term, new_term->minimal_polynomial);
        increment_reference_count(new_term);
        struct Number*out = sum_incorporate_term_in_terms_of_new_generator(pool_set, stack_a,
            stack_b, a, new_term, new_term, &a_terms_in_terms_of_new_term, x);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    struct Number*new_generator = number_allocate(pool_set);
    new_generator->operation = 'g';
    new_generator->term_count = a->generator->term_count + 1;
    new_generator->terms =
        pool_value_allocate(pool_set, new_generator->term_count * sizeof(struct Number*));
    memcpy(new_generator->terms, a->generator->terms,
        a->generator->term_count * sizeof(struct Number*));
    struct Rational*k = &rational_one;
    while (true)
    {
        increment_reference_count(new_term);
        new_generator->terms[a->generator->term_count] = number_rational_multiply(pool_set, stack_a,
            stack_b, new_term, k);
        struct RationalPolynomial*term_in_terms_of_new_generator =
            number_a_in_terms_of_b(pool_set, stack_a, stack_b, new_term, new_generator);
        if (term_in_terms_of_new_generator)
        {
            number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, new_generator);
            struct Matrix a_terms_in_terms_of_new_generator;
            sun_convert_terms_to_new_generator(stack_a, stack_b, &a_terms_in_terms_of_new_generator,
                a, rational_polynomial_subtract(stack_a, stack_b, x,
                    rational_polynomial_rational_multiply(stack_a, stack_b,
                        term_in_terms_of_new_generator, k)), new_generator->minimal_polynomial);
            struct Number*out = sum_incorporate_term_in_terms_of_new_generator(pool_set, stack_a,
                stack_b, a, new_term, new_generator, &a_terms_in_terms_of_new_generator,
                term_in_terms_of_new_generator);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        number_free(pool_set, new_generator->terms[a->generator->term_count]);
        k->numerator = integer_add(stack_a, k->numerator, &one);
    }
}

struct Number*number_add(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a, struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->value_count == 0)
        {
            number_free(pool_set, a);
            return b;
        }
        if (b->operation == 'r')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct Number*out = number_allocate(pool_set);
            out->operation = 'r';
            out->value = *rational_add(stack_a, stack_b, &a->value, &b->value);
            rational_move_value_to_pool(pool_set, &out->value);
            number_free(pool_set, a);
            number_free(pool_set, b);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        break;
    case '^':
    case '*':
        switch (b->operation)
        {
        case 'r':
        {
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            increment_reference_count(a);
            out->generator = a;
            out->first_term = pool_value_allocate(pool_set, sizeof(struct Term));
            out->first_term->value = a;
            out->first_term->in_terms_of_generator = pool_polynomial_allocate(pool_set, 2);
            out->first_term->in_terms_of_generator->coefficients[0] =
                pool_value_allocate(pool_set, sizeof(struct Rational));
            out->first_term->in_terms_of_generator->coefficients[0] =
                rational_copy_to_pool(pool_set, &rational_zero);
            out->first_term->in_terms_of_generator->coefficients[1] =
                rational_copy_to_pool(pool_set, &rational_one);
            out->first_term->next = pool_value_allocate(pool_set, sizeof(struct Term));
            out->first_term->next->value = b;
            out->first_term->next->in_terms_of_generator = pool_polynomial_allocate(pool_set, 1);
            out->first_term->next->in_terms_of_generator->coefficients[0] =
                rational_copy_to_pool(pool_set, &b->value);
            out->first_term->next->next = 0;
            return out;
        }
        case '^':
        case '*':
        {
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            increment_reference_count(a);
            out->generator = a;
            out->first_term = pool_value_allocate(pool_set, sizeof(struct Term));
            out->first_term->value = a;
            out->first_term->in_terms_of_generator = pool_polynomial_allocate(pool_set, 2);
            out->first_term->in_terms_of_generator->coefficients[0] =
                pool_value_allocate(pool_set, sizeof(struct Rational));
            out->first_term->in_terms_of_generator->coefficients[0] =
                rational_copy_to_pool(pool_set, &rational_zero);
            out->first_term->in_terms_of_generator->coefficients[1] =
                rational_copy_to_pool(pool_set, &rational_one);
            out->first_term->next = 0;
            return sum_incorporate_term(pool_set, stack_a, stack_b, out, b);
        }
        }
        break;
    case '+':
        if (b->operation == '+')
        {
            struct Term*term = b->first_term;
            while (term)
            {
                increment_reference_count(term->value);
                a = number_add(pool_set, stack_a, stack_b, a, term->value);
                term = term->next;
            }
            number_free(pool_set, b);
            return a;
        }
        else
        {
            return sum_incorporate_term(pool_set, stack_a, stack_b, a, b);
        }
    }
    return number_add(pool_set, stack_a, stack_b, b, a);
}

struct Number*number_rational_multiply(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*b)
{
    if (b->numerator->value_count == 0)
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value.numerator = &zero;
        out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
        number_free(pool_set, a);
        return out;
    }
    if (integer_equals(b->numerator, &one) && integer_equals(b->denominator, &one))
    {
        return a;
    }
    switch (a->operation)
    {
    case 'r':
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value = *rational_multiply(stack_a, stack_b, &a->value, b);
        rational_move_value_to_pool(pool_set, &out->value);
        number_free(pool_set, a);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '^':
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = '*';
        out->left = number_allocate(pool_set);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy_to_pool(pool_set, b->numerator);
        out->left->value.denominator = integer_copy_to_pool(pool_set, b->denominator);
        out->right = a;
        return out;
    }
    case '*':
    {
        if (a->left->operation == 'r')
        {
            increment_reference_count(a->left);
            increment_reference_count(a->right);
            struct Number*out = number_multiply(pool_set, stack_a, stack_b,
                number_rational_multiply(pool_set, stack_a, stack_b, a->left, b), a->right);
            number_free(pool_set, a);
            return out;
        }
        struct Number*out = number_allocate(pool_set);
        out->operation = '*';
        out->left = number_allocate(pool_set);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy_to_pool(pool_set, b->numerator);
        out->left->value.denominator = integer_copy_to_pool(pool_set, b->denominator);
        out->right = a;
        return out;
    }
    case '+':
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct Number*out = number_allocate(pool_set);
        out->operation = '+';
        increment_reference_count(a->generator);
        out->generator = a->generator;
        struct Term**a_term = &a->first_term;
        struct Term**out_term = &out->first_term;
        while (*a_term)
        {
            increment_reference_count(*a_term);
            (*out_term)->value =
                number_rational_multiply(pool_set, stack_a, stack_b, (*a_term)->value, b);
            (*out_term)->in_terms_of_generator = rational_polynomial_copy_to_pool(pool_set,
                rational_polynomial_rational_multiply(stack_a, stack_b,
                    (*a_term)->in_terms_of_generator, b));
            stack_a->cursor = stack_a_savepoint;
            *a_term = (*a_term)->next;
            *out_term = (*out_term)->next;
        }
        number_free(pool_set, a);
        return out;
    }
    case 'g':
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = 'g';
        out->term_count = a->term_count;
        out->terms = pool_value_allocate(pool_set, out->term_count * sizeof(struct Number*));
        for (size_t i = 0; i < out->term_count; ++i)
        {
            increment_reference_count(a->terms[i]);
            out->terms[i] = number_rational_multiply(pool_set, stack_a, stack_b, a->terms[i], b);
        }
        number_free(pool_set, a);
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_consolidate_product(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Number*b)
{
    if (b->operation == 'r')
    {
        struct Number*out = number_rational_multiply(pool_set, stack_a, stack_b, a, &b->value);
        number_free(pool_set, b);
        return out;
    }
    switch (a->operation)
    {
    case 'r':
    {
        struct Number*out = number_rational_multiply(pool_set, stack_a, stack_b, b, &a->value);
        number_free(pool_set, a);
        return out;
    }
    case '^':
    {
        if (b->operation == '^')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct ExtendedGCDInfo info;
            integer_extended_gcd(stack_a, stack_b, &info, a->right->value.denominator,
                b->right->value.denominator);
            struct Integer*product_index = integer_euclidean_quotient(stack_a, stack_b,
                integer_multiply(stack_a, stack_b, a->right->value.denominator,
                    b->right->value.denominator), info.gcd);
            increment_reference_count(a->left);
            struct Number*radicand_factor_a = number_exponentiate(pool_set, stack_a, stack_b,
                a->left, &(struct Rational){info.b_over_gcd, &one});
            increment_reference_count(b->left);
            struct Number*radicand_factor_b = number_exponentiate(pool_set, stack_a, stack_b,
                b->left, &(struct Rational){info.a_over_gcd, &one});
            struct RationalInterval radicand_factor_a_argument_estimate;
            number_rational_argument_estimate(pool_set, stack_a, stack_b,
                &radicand_factor_a_argument_estimate, radicand_factor_a, &rational_one);
            struct RationalInterval radicand_factor_b_argument_estimate;
            number_rational_argument_estimate(pool_set, stack_a, stack_b,
                &radicand_factor_b_argument_estimate, radicand_factor_b, &rational_one);
            if (!integer_equals(a->right->value.denominator, b->right->value.denominator))
            {
                struct Rational*radicand_factor_argument_interval_size =
                    rational_integer_divide(stack_a, stack_b, &pi_estimate_min, product_index);
                struct RationalInterval a_argument_estimate;
                number_rational_argument_estimate(pool_set, stack_a, stack_b, &a_argument_estimate,
                    a, radicand_factor_argument_interval_size);
                struct RationalInterval b_argument_estimate;
                number_rational_argument_estimate(pool_set, stack_a, stack_b, &b_argument_estimate,
                    b, radicand_factor_argument_interval_size);
                if (rational_compare(stack_a, stack_b,
                    rational_integer_divide(stack_a, stack_b,
                        radicand_factor_a_argument_estimate.max, product_index),
                    a_argument_estimate.min) < 0 ||
                    rational_compare(stack_a, stack_b,
                        rational_integer_divide(stack_a, stack_b,
                            radicand_factor_b_argument_estimate.max, product_index),
                        b_argument_estimate.min) < 0)
                {
                    number_free(pool_set, radicand_factor_a);
                    number_free(pool_set, radicand_factor_b);
                    stack_a->cursor = stack_a_savepoint;
                    return number_product_consolidation_failed;
                }
            }
            struct Number*product_radicand =
                number_multiply(pool_set, stack_a, stack_b, radicand_factor_a, radicand_factor_b);
            struct Number*out = number_exponentiate(pool_set, stack_a, stack_b, product_radicand,
                &(struct Rational){&one, product_index});
            struct RationalInterval product_radicand_rational_argument_estimate;
            number_rational_argument_estimate(pool_set, stack_a, stack_b,
                &product_radicand_rational_argument_estimate, product_radicand, &pi_estimate_min);
            if (rational_compare(stack_a, stack_b, product_radicand_rational_argument_estimate.max,
                rational_add(stack_a, stack_b, radicand_factor_a_argument_estimate.min,
                    radicand_factor_b_argument_estimate.min)) < 0)
            {
                struct Number*root_of_unity =
                    get_roots_of_unity(stack_a, stack_b, product_index)[1];
                increment_reference_count(root_of_unity);
                out = number_multiply(pool_set, stack_a, stack_b, out, root_of_unity);
            }
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        break;
    }
    case '*':
        switch (b->operation)
        {
        case '^':
        {
            struct Number*factor =
                number_consolidate_product(pool_set, stack_a, stack_b, a->left, b);
            if (factor)
            {
                return number_multiply(pool_set, stack_a, stack_b, a->right, factor);
            }
            factor = number_consolidate_product(pool_set, stack_a, stack_b, a->right, b);
            if (factor)
            {
                return number_multiply(pool_set, stack_a, stack_b, a->left, factor);
            }
            return number_product_consolidation_failed;
        }
        case '*':
        {
            struct Number*out = number_multiply(pool_set, stack_a, stack_b, a->left,
                number_multiply(pool_set, stack_a, stack_b, a->right, b));
            a->left = 0;
            a->right = 0;
            number_free(pool_set, a);
            return out;
        }
        }
        break;
    case '+':
    case 'g':
    {
        struct TermIterator iterator;
        term_iterator_initialize(&iterator, a);
        increment_reference_count(iterator.term);
        increment_reference_count(b);
        struct Number*out = number_multiply(pool_set, stack_a, stack_b, iterator.term, b);
        term_iterator_increment(&iterator);
        while (iterator.term)
        {
            increment_reference_count(iterator.term);
            increment_reference_count(b);
            out = number_add(pool_set, stack_a, stack_b, out,
                number_multiply(pool_set, stack_a, stack_b, iterator.term, b));
            term_iterator_increment(&iterator);
        }
        number_free(pool_set, a);
        number_free(pool_set, b);
        return out;
    }
    }
    return number_multiply(pool_set, stack_a, stack_b, b, a);
}

struct Number*number_multiply(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a, struct Number*b)
{
    struct Number*out = number_consolidate_product(pool_set, stack_a, stack_b, a, b);
    if (out == number_product_consolidation_failed)
    {
        out = number_allocate(pool_set);
        out->operation = '*';
        out->left = a;
        out->right = b;
    }
    return out;
}

struct Number*number_reciprocal(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        if (!a->value.numerator->value_count)
        {
            puts("Tried to divide by 0.");
            return number_divide_by_zero_error;
        }
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        increment_reference_count(a->value.numerator);
        increment_reference_count(a->value.denominator);
        out->value.numerator = a->value.denominator;
        out->value.denominator = a->value.numerator;
        number_free(pool_set, a);
        return out;
    }
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        increment_reference_count(a->left);
        increment_reference_count(a->left);
        struct Number*out = number_divide(pool_set, stack_a, stack_b,
            number_exponentiate(pool_set, stack_a, stack_b, a->left,
                &(struct Rational){ integer_add(stack_a, a->right->value.denominator, &INT(1, -)),
                a->right->value.denominator }), a->left);
        number_free(pool_set, a);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '*':
    {
        increment_reference_count(a->left);
        increment_reference_count(a->right);
        struct Number*out = number_multiply(pool_set, stack_a, stack_b,
            number_reciprocal(pool_set, stack_a, stack_b, a->left),
            number_reciprocal(pool_set, stack_a, stack_b, a->right));
        number_free(pool_set, a);
        return out;
    }
    case '+':
    {
        void*stack_a_savepoint = stack_a->cursor;
        number_calculate_conjugates(pool_set, stack_a, stack_b, a);
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value.numerator = pool_integer_initialize(pool_set, 1, 1);
        increment_reference_count(out->value.numerator);
        out->value.denominator = out->value.numerator;
        struct Number*conjugate = a->next;
        while (conjugate)
        {
            increment_reference_count(conjugate);
            out = number_multiply(pool_set, stack_a, stack_b, out, conjugate);
            conjugate = conjugate->next;
        }
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a);
        struct Rational*coefficient =
            rational_reciprocal(stack_a, stack_b, a->minimal_polynomial->coefficients[0]);
        if (a->minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient->numerator->sign *= -1;
        }
        out = number_rational_multiply(pool_set, stack_a, stack_b, out, coefficient);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_divide(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*dividend, struct Number*divisor)
{
    struct Number*reciprocal = number_reciprocal(pool_set, stack_a, stack_b, divisor);
    if (reciprocal == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    return number_multiply(pool_set, stack_a, stack_b, dividend, reciprocal);
}

struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        return &a->value;
    case '^':
        return &rational_one;
    case '*':
        return number_rational_factor(output_stack, local_stack, a->left);
    case '+':
    case 'g':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct TermIterator iterator;
        term_iterator_initialize(&iterator, a);
        struct Rational*out = &rational_one;
        while (iterator.term)
        {
            out = rational_reduced(local_stack, output_stack,
                integer_gcd(local_stack, output_stack, out->numerator,
                    number_rational_factor(local_stack, output_stack, iterator.term)->numerator),
                integer_lcm(local_stack, output_stack, out->denominator,
                    number_rational_factor(local_stack, output_stack, iterator.term)->denominator));
            term_iterator_increment(&iterator);
        }
        out = rational_copy_to_stack(output_stack, out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.\n");
    }
}

struct Number*number_exponentiate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*base, struct Rational*exponent)
{
    if (exponent->numerator->sign < 0)
    {
        base = number_reciprocal(pool_set, stack_a, stack_b, base);
        if (base == number_divide_by_zero_error)
        {
            return number_divide_by_zero_error;
        }
        void*stack_a_savepoint = stack_a->cursor;
        struct Number*out = number_exponentiate(pool_set, stack_a, stack_b, base,
            rational_negative(stack_a, exponent));
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    switch (base->operation)
    {
    case 'r':
    {
        if (base->value.numerator->value_count == 0)
        {
            struct Number*out = number_allocate(pool_set);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
            return out;
        }
        void*stack_a_savepoint = stack_a->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*out = number_allocate(pool_set);
            out->operation = 'r';
            out->value =
                *rational_exponentiate(stack_a, stack_b, &base->value, exponent->numerator);
            rational_move_value_to_pool(pool_set, &out->value);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        if (!integer_equals(base->value.denominator, &one))
        {
            struct Number*new_numerator_base = number_allocate(pool_set);
            new_numerator_base->operation = 'r';
            new_numerator_base->value.denominator = pool_integer_initialize(pool_set, 1, 1);
            new_numerator_base->value.numerator =
                integer_multiply(stack_a, stack_b, base->value.numerator,
                    integer_exponentiate(stack_a, stack_b, base->value.denominator,
                        integer_add(stack_a, exponent->denominator, &INT(1, -))));
            integer_move_to_pool(pool_set, &new_numerator_base->value.numerator);
            struct Number*new_numerator =
                number_exponentiate(pool_set, stack_a, stack_b, new_numerator_base, exponent);
            struct Number*out = number_rational_multiply(pool_set, stack_a, stack_b, new_numerator,
                &(struct Rational){&one, base->value.denominator});
            number_free(pool_set, base);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        if (!integer_equals(exponent->denominator, &one))
        {
            struct Number*radicand = number_allocate(pool_set);
            radicand->operation = 'r';
            radicand->value.numerator = integer_magnitude(stack_a, base->value.numerator);
            radicand->value.denominator = pool_integer_initialize(pool_set, 1, 1);
            struct Factor*factors;
            size_t factor_count =
                integer_factor(stack_a, stack_b, &factors, radicand->value.numerator);
            struct Integer*coefficient = &one;
            struct Integer*multiplicity_gcd = exponent->denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct IntegerDivision multiplicity_reduction;
                integer_euclidean_divide(stack_a, stack_b, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(stack_a, stack_b, coefficient,
                    integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = integer_gcd(stack_a, stack_b, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            if (base->value.numerator->sign > 0)
            {
                struct Integer*reduced_degree = integer_euclidean_quotient(stack_a, stack_b,
                    exponent->denominator, multiplicity_gcd);
                if (integer_equals(reduced_degree, &one))
                {
                    radicand->value.numerator = integer_copy_to_pool(pool_set, coefficient);
                    stack_a->cursor = stack_a_savepoint;
                    return radicand;
                }
                exponent = STACK_SLOT_ALLOCATE(stack_a, struct Rational);
                exponent->numerator = &one;
                exponent->denominator = reduced_degree;
            }
            radicand->value.numerator =
                stack_integer_initialize(stack_a, 1, base->value.numerator->sign);
            number_free(pool_set, base);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Integer*reduced_multiplicity = integer_euclidean_quotient(stack_a, stack_b,
                    factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(stack_a, stack_b,
                    factors[factor_index].value, reduced_multiplicity);
                radicand->value.numerator =
                    integer_multiply(stack_a, stack_b, radicand->value.numerator, exponentiation);
            }
            integer_move_to_pool(pool_set, &radicand->value.numerator);
            struct Number*number_exponent = number_allocate(pool_set);
            number_exponent->operation = 'r';
            number_exponent->value.numerator = pool_integer_initialize(pool_set, 1, 1);
            number_exponent->value.denominator =
                integer_copy_to_pool(pool_set, exponent->denominator);
            struct Number*surd = number_allocate(pool_set);
            surd->operation = '^';
            surd->left = radicand;
            surd->right = number_exponent;
            stack_a->cursor = stack_a_savepoint;
            return number_rational_multiply(pool_set, stack_a, stack_b, surd,
                &(struct Rational){coefficient, &one});
        }
        return base;
    }
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        exponent = rational_multiply(stack_a, stack_b, exponent, &base->right->value);
        increment_reference_count(base->left);
        struct Number*out = number_exponentiate(pool_set, stack_a, stack_b, base->left, exponent);
        number_free(pool_set, base);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '*':
    {
        increment_reference_count(base->left);
        increment_reference_count(base->right);
        struct Number*out = number_multiply(pool_set, stack_a, stack_b,
            number_exponentiate(pool_set, stack_a, stack_b, base->left, exponent),
            number_exponentiate(pool_set, stack_a, stack_b, base->right, exponent));
        number_free(pool_set, base);
        return out;
    }
    case '+':
    {
        void*stack_a_savepoint = stack_a->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*exponentiation_by_numerator = number_allocate(pool_set);
            exponentiation_by_numerator->operation = 'r';
            exponentiation_by_numerator->value.numerator = pool_integer_initialize(pool_set, 1, 1);
            exponentiation_by_numerator->value.denominator =
                pool_integer_initialize(pool_set, 1, 1);
            struct Integer*numerator = exponent->numerator;
            while (numerator->sign > 0)
            {
                if (numerator->value[0] & 1)
                {
                    increment_reference_count(base);
                    exponentiation_by_numerator = number_multiply(pool_set, stack_a, stack_b,
                        exponentiation_by_numerator, base);
                }
                increment_reference_count(base);
                base = number_multiply(pool_set, stack_a, stack_b, base, base);
                numerator = integer_half(stack_a, numerator);
            }
            number_free(pool_set, base);
            stack_a->cursor = stack_a_savepoint;
            return number_exponentiate(pool_set, stack_a, stack_b, exponentiation_by_numerator,
                &(struct Rational){&one, exponent->denominator});
        }
        struct Rational*base_rational_factor = number_rational_factor(stack_a, stack_b, base);
        struct Factor*factors;
        size_t factor_count =
            integer_factor(stack_a, stack_b, &factors, base_rational_factor->numerator);
        struct Rational base_cancelling_rational_factor =
            { integer_exponentiate(stack_a, stack_b, base_rational_factor->denominator,
                exponent->denominator), &one };
        struct Rational product_rational_factor = { &one, base_rational_factor->denominator };
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_quotient(stack_a, stack_b,
                factors[factor_index].multiplicity, exponent->denominator);
            base_cancelling_rational_factor.denominator = integer_multiply(stack_a, stack_b,
                integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                    integer_multiply(stack_a, stack_b, reduced_multiplicity,
                        exponent->denominator)), base_cancelling_rational_factor.denominator);
            product_rational_factor.numerator =
                integer_multiply(stack_a, stack_b, product_rational_factor.numerator,
                    integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                        reduced_multiplicity));
        }
        struct Number*surd = number_allocate(pool_set);
        surd->operation = '^';
        surd->left = number_rational_multiply(pool_set, stack_a, stack_b, base,
            &base_cancelling_rational_factor);
        surd->right = number_allocate(pool_set);
        surd->right->operation = 'r';
        surd->right->value.numerator = pool_integer_initialize(pool_set, 1, 1);
        surd->right->value.denominator = integer_copy_to_pool(pool_set, exponent->denominator);
        struct Number*out =
            number_rational_multiply(pool_set, stack_a, stack_b, surd, &product_rational_factor);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}