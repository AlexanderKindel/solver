#include "declarations.h"

struct Number*sum_incorporate_term_in_terms_of_new_generator(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Number*new_term,
    struct Number*new_generator, struct Matrix*a_terms_in_terms_of_new_generator,
    struct RationalPolynomial*new_term_in_terms_of_new_generator)
{
    if (new_term_in_terms_of_new_generator->coefficient_count >
        a_terms_in_terms_of_new_generator->height)
    {
        struct Number*out = ALLOCATE(output_stack, struct Number);
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
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out = ALLOCATE(output_stack, struct Number);
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
            out->terms_in_terms_of_generator[a->term_count] = new_term_in_terms_of_new_generator;
            out = number_copy(output_stack, out);
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
        struct Number*out = ALLOCATE(output_stack, struct Number);
        out->operation = 'r';
        out->value = rational_zero;
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case 1:
    {
        out = number_copy(output_stack, out->terms[0]);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        out = number_copy(output_stack, out);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
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

struct Number*sum_incorporate_term(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Number*new_term)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalPolynomial*a_generator_minimal_polynomial =
        number_minimal_polynomial(local_stack, output_stack, a->generator);
    struct RationalPolynomial*new_term_minimal_polynomial =
        number_minimal_polynomial(local_stack, output_stack, new_term);
    struct RationalPolynomial*new_term_in_terms_of_a_generator =
        number_a_in_terms_of_b(local_stack, output_stack, new_term, new_term_minimal_polynomial,
            a->generator, a_generator_minimal_polynomial);
    if (new_term_in_terms_of_a_generator)
    {
        struct Matrix a_terms_in_terms_of_a_generator = { ARRAY_ALLOCATE(local_stack,
            a_generator_minimal_polynomial->coefficient_count - 1, struct Rational**),
            a->term_count, a_generator_minimal_polynomial->coefficient_count - 1 };
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
        number_a_in_terms_of_b(local_stack, output_stack, a->generator,
            a_generator_minimal_polynomial, new_term, new_term_minimal_polynomial);
    if (old_generator_in_terms_of_new_term)
    {
        struct RationalPolynomial*new_term_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, new_term);
        struct Matrix a_terms_in_terms_of_new_term;
        sun_convert_terms_to_new_generator(local_stack, output_stack, &a_terms_in_terms_of_new_term,
            a, old_generator_in_terms_of_new_term, new_term_minimal_polynomial);
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
        struct RationalPolynomial*new_generator_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, new_generator);
        struct RationalPolynomial*term_in_terms_of_new_generator =
            number_a_in_terms_of_b(local_stack, output_stack, new_term, new_term_minimal_polynomial,
                new_generator, new_generator_minimal_polynomial);
        if (term_in_terms_of_new_generator)
        {
            struct Matrix a_terms_in_terms_of_new_generator;
            sun_convert_terms_to_new_generator(local_stack, output_stack,
                &a_terms_in_terms_of_new_generator, a,
                rational_polynomial_subtract(local_stack, output_stack,
                    (struct RationalPolynomial*)&x,
                    rational_polynomial_rational_multiply(local_stack, output_stack,
                        term_in_terms_of_new_generator, k)), new_generator_minimal_polynomial);
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
            struct Number*out = ALLOCATE(output_stack, struct Number);
            out->operation = 'r';
            out->value = *rational_add(local_stack, output_stack, &a->value, &b->value);
            out->value.numerator = integer_copy(output_stack, out->value.numerator);
            out->value.denominator = integer_copy(output_stack, out->value.denominator);
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

struct Number*number_rational_multiply(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Rational*b)
{
    if (b->numerator->value_count == 0)
    {
        struct Number*out = ALLOCATE(output_stack, struct Number);
        out->operation = 'r';
        out->value = rational_zero;
        return out;
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
        out = ALLOCATE(output_stack, struct Number);
        out->operation = 'r';
        out->value = *rational_multiply(local_stack, output_stack, &a->value, b);
        out->value.numerator = integer_copy(output_stack, out->value.numerator);
        out->value.denominator = integer_copy(output_stack, out->value.denominator);
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
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '*';
        out->left = ALLOCATE(output_stack, struct Number);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy(output_stack, b->numerator);
        out->left->value.denominator = integer_copy(output_stack, b->denominator);
        out->right = number_copy(output_stack, a);
        return out;
    case '+':
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '+';
        out->generator = number_copy(output_stack, a->generator);
        out->terms_in_terms_of_generator =
            ARRAY_ALLOCATE(output_stack, a->term_count, struct RationalPolynomial*);
        for (size_t i = 0; i < a->term_count; ++i)
        {
            out->terms_in_terms_of_generator[i] =
                rational_polynomial_rational_multiply(output_stack, local_stack,
                    a->terms_in_terms_of_generator[i], b);
        }
        break;
    case 'g':
    {
        out = ALLOCATE(output_stack, struct Number);
        out->operation = 'g';
        break;
    }
    default:
        crash("Number operation not recognized.");
    }
    out->term_count = a->term_count;
    out->terms = ARRAY_ALLOCATE(output_stack, a->term_count, struct Number*);
    for (size_t i = 0; i < a->term_count; ++i)
    {
        out->terms[i] = number_rational_multiply(output_stack, local_stack, a->terms[i], b);
    }
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
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '*';
        out->left = number_copy(output_stack, a);
        out->right = number_copy(output_stack, b);
    }
    return out;
}

struct Number*number_reciprocal(struct Stack*output_stack, struct Stack*local_stack,
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
        struct Number*out = ALLOCATE(output_stack, struct Number);
        out->operation = 'r';
        out->value.numerator = integer_copy(output_stack, a->value.denominator);
        out->value.denominator = integer_copy(output_stack, a->value.numerator);
        return out;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_divide(output_stack, local_stack,
            number_exponentiate(local_stack, output_stack, a->left,
                &(struct Rational){ integer_add(local_stack, a->right->value.denominator,
                    &INT(1, -)), a->right->value.denominator }), a->left);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_multiply(output_stack, local_stack,
            number_reciprocal(local_stack, output_stack, a->left),
            number_reciprocal(local_stack, output_stack, a->right));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalPolynomial*minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a);
        struct Number**conjugates =
            number_conjugates(local_stack, output_stack, a, minimal_polynomial);
        struct Number*out = &number_one;
        for (size_t i = 0; i < minimal_polynomial->coefficient_count - 1; ++i)
        {
            out = number_multiply(local_stack, output_stack, out, conjugates[i]);
        }
        struct Rational*coefficient =
            rational_reciprocal(local_stack, output_stack, minimal_polynomial->coefficients[0]);
        if (minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient->numerator->sign *= -1;
        }
        out = number_rational_multiply(output_stack, local_stack, out, coefficient);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
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
            struct Number*out = ALLOCATE(output_stack, struct Number);
            out->operation = 'r';
            out->value = *rational_exponentiate(local_stack, output_stack, &base->value,
                exponent->numerator);
            out->value.numerator = integer_copy(output_stack, out->value.numerator);
            out->value.denominator = integer_copy(output_stack, out->value.denominator);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (!integer_equals(base->value.denominator, &one))
        {
            struct Number*new_numerator_base = ALLOCATE(local_stack, struct Number);
            new_numerator_base->operation = 'r';
            new_numerator_base->value.denominator = integer_initialize(local_stack, 1, 1);
            new_numerator_base->value.numerator =
                integer_multiply(local_stack, output_stack, base->value.numerator,
                    integer_exponentiate(local_stack, output_stack, base->value.denominator,
                        integer_add(local_stack, exponent->denominator, &INT(1, -))));
            struct Number*out = number_rational_multiply(output_stack, local_stack,
                number_exponentiate(local_stack, output_stack, new_numerator_base, exponent),
                &(struct Rational){&one, base->value.denominator});
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        if (!integer_equals(exponent->denominator, &one))
        {
            struct Number*radicand = ALLOCATE(local_stack, struct Number);
            radicand->operation = 'r';
            radicand->value = *rational_magnitude(local_stack, &base->value);
            struct Factor*factors;
            size_t factor_count =
                integer_factor(local_stack, output_stack, &factors, radicand->value.numerator);
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
            struct Number*number_exponent = ALLOCATE(local_stack, struct Number);
            number_exponent->operation = 'r';
            if (base->value.numerator->sign > 0)
            {
                struct Integer*reduced_degree = integer_euclidean_quotient(local_stack,
                    output_stack, exponent->denominator, multiplicity_gcd);
                if (integer_equals(reduced_degree, &one))
                {
                    radicand->value.numerator = coefficient;
                    local_stack->cursor = local_stack_savepoint;
                    return number_copy(output_stack, radicand);
                }
                number_exponent->value.numerator = &one;
                number_exponent->value.denominator = reduced_degree;
            }
            else
            {
                number_exponent->value = *exponent;
            }
            radicand->value.numerator =
                integer_initialize(local_stack, 1, base->value.numerator->sign);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Integer*reduced_multiplicity = integer_euclidean_quotient(local_stack,
                    output_stack, factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(local_stack, output_stack,
                    factors[factor_index].value, reduced_multiplicity);
                radicand->value.numerator = integer_multiply(local_stack, output_stack,
                    radicand->value.numerator, exponentiation);
            }
            struct Number*surd = ALLOCATE(local_stack, struct Number);
            surd->operation = '^';
            surd->left = radicand;
            surd->right = number_exponent;
            struct Number*out = number_rational_multiply(output_stack, local_stack, surd,
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
        struct Number*surd = ALLOCATE(local_stack, struct Number);
        surd->operation = '^';
        surd->left = number_rational_multiply(local_stack, output_stack, base,
            &base_cancelling_rational_factor);
        surd->right = ALLOCATE(local_stack, struct Number);
        surd->right->operation = 'r';
        surd->right->value = *exponent;
        struct Number*out =
            number_rational_multiply(output_stack, local_stack, surd, &product_rational_factor);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}