#include "declarations.h"

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

struct RationalPolynomial*number_minimal_polynomial_from_algebraic_form(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct AlgebraicNumber*algebraic_form,
    struct RationalPolynomial*generator_minimal_polynomials[2])
{
    void*local_stack_savepoint = local_stack->cursor;
    size_t term_count_upper_bound = (generator_minimal_polynomials[0]->coefficient_count - 1) *
        (generator_minimal_polynomials[1]->coefficient_count - 1) + 1;
    struct Matrix matrix;
    matrix.rows = ARRAY_ALLOCATE(local_stack, term_count_upper_bound, struct Rational**);
    matrix.height = 0;
    size_t(**terms_present)[2] = ARRAY_ALLOCATE(local_stack, term_count_upper_bound, size_t(*)[2]);
    size_t terms_present_count = 1;
    struct AlgebraicNumber*power = ALLOCATE(local_stack, struct AlgebraicNumber);
    power->term_coefficient = &rational_one;
    power->generator_degrees[0] = 0;
    power->generator_degrees[1] = 0;
    power->next_term = 0;
    struct AlgebraicNumber**powers =
        ARRAY_ALLOCATE(local_stack, term_count_upper_bound, struct AlgebraicNumber*);
    bool constant_is_present = false;
    while (!constant_is_present || matrix.height < terms_present_count)
    {
        power = algebraic_number_multiply(local_stack, output_stack, power, algebraic_form,
            generator_minimal_polynomials);
        powers[matrix.height] = power;
        struct AlgebraicNumber*power_term = powers[matrix.height];
        matrix.rows[matrix.height] =
            ARRAY_ALLOCATE(local_stack, term_count_upper_bound, struct Rational*);
        for (size_t i = 0; i < term_count_upper_bound; ++i)
        {
            matrix.rows[matrix.height][i] = &rational_zero;
        }
        while (power_term)
        {
            for (size_t i = 0; i < terms_present_count; ++i)
            {
                if (degree_compare(*terms_present[i], power_term->generator_degrees) == 0)
                {
                    matrix.rows[matrix.height][i] = power_term->term_coefficient;
                    goto degree_already_present;
                }
            }
            terms_present[terms_present_count] = &power_term->generator_degrees;
            matrix.rows[matrix.height][matrix.height] = power_term->term_coefficient;
        degree_already_present:
            power_term = power_term->next_term;
        }
        if (matrix.rows[matrix.height][0]->numerator->value_count)
        {
            constant_is_present = true;
        }
        ++matrix.height;
    }
    matrix.width = matrix.height;
    struct RationalPolynomial**augmentation =
        ARRAY_ALLOCATE(local_stack, matrix.height, struct RationalPolynomial*);
    for (size_t i = 0; i < matrix.height; ++i)
    {
        array_reverse(matrix.rows[i], matrix.width);
        augmentation[i] = polynomial_allocate(local_stack, i + 2);
        for (size_t j = 0; j <= i; ++j)
        {
            augmentation[i]->coefficients[j] = &rational_zero;
        }
        augmentation[i]->coefficients[i + 1] = &rational_one;
    }
    matrix_row_echelon_form(rational_polynomial_rational_multiply, rational_polynomial_subtract,
        local_stack, output_stack, &matrix, augmentation);
    struct RationalPolynomial*annulling_polynomial = augmentation[matrix.height - 1];
    annulling_polynomial->coefficients[0] = rational_subtract(local_stack, output_stack,
        annulling_polynomial->coefficients[0], matrix.rows[matrix.height - 1][matrix.width - 1]);
    struct RationalPolynomial*out =
        number_minimal_polynomial_from_annulling_polynomial(output_stack, local_stack,
            annulling_polynomial, a);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct RationalPolynomial*number_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        struct RationalPolynomial*out = polynomial_allocate(output_stack, 2);
        out->coefficients[0] = rational_negative(output_stack, &a->value);
        out->coefficients[1] = &rational_one;
        return out;
    }
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalPolynomial*radicand_minimal_polynomial =
            number_minimal_polynomial(output_stack, local_stack, a->left);
        size_t surd_index = integer_to_size_t(a->left->value.denominator);
        struct RationalPolynomial*annulling_polynomial = polynomial_allocate(local_stack,
            surd_index * (radicand_minimal_polynomial->coefficient_count - 1) + 1);
        annulling_polynomial->coefficients[0] = radicand_minimal_polynomial->coefficients[0];
        for (size_t i = 1; i < radicand_minimal_polynomial->coefficient_count; ++i)
        {
            size_t annulling_polynomial_coefficient_index = i * surd_index;
            annulling_polynomial->coefficients[annulling_polynomial_coefficient_index] =
                radicand_minimal_polynomial->coefficients[i];
            for (size_t j = 1; j < surd_index; ++j)
            {
                annulling_polynomial->coefficients[annulling_polynomial_coefficient_index - j] =
                    &rational_zero;
            }
        }
        struct RationalPolynomial*out =
            number_minimal_polynomial_from_annulling_polynomial(output_stack, local_stack,
                annulling_polynomial, a);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct AlgebraicNumber algebraic_form;
        algebraic_form.term_coefficient = &rational_one;
        algebraic_form.generator_degrees[0] = 1;
        algebraic_form.generator_degrees[1] = 1;
        algebraic_form.next_term = 0;
        struct RationalPolynomial*out = number_minimal_polynomial_from_algebraic_form(output_stack,
            local_stack, a, &algebraic_form, (struct RationalPolynomial*[2]) {
                number_minimal_polynomial(local_stack, output_stack, a->left),
                number_minimal_polynomial(local_stack, output_stack, a->right) });
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        if (a->term_count == 1)
        {
            return number_minimal_polynomial(output_stack, local_stack, a->terms[0]);
        }
        void*local_stack_savepoint = local_stack->cursor;
        --a->term_count;
        struct RationalPolynomial*left_terms_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a);
        ++a->term_count;
        struct AlgebraicNumber left_terms_algebraic_form;
        left_terms_algebraic_form.term_coefficient = &rational_one;
        left_terms_algebraic_form.generator_degrees[0] = 1;
        left_terms_algebraic_form.generator_degrees[1] = 0;
        left_terms_algebraic_form.next_term = 0;
        struct AlgebraicNumber algebraic_form;
        algebraic_form.term_coefficient = &rational_one;
        algebraic_form.generator_degrees[0] = 0;
        algebraic_form.generator_degrees[1] = 1;
        algebraic_form.next_term = &left_terms_algebraic_form;
        struct RationalPolynomial*out =
            number_minimal_polynomial_from_algebraic_form(output_stack, local_stack, a,
                &algebraic_form, (struct RationalPolynomial*[2]){left_terms_minimal_polynomial,
                number_minimal_polynomial(local_stack, output_stack, a->terms[a->term_count - 1])});       
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
            struct RationalPolynomial*candidate_minimal_polynomial =
                number_minimal_polynomial(local_stack, output_stack, candidate);
            if (rational_polynomial_equals(minimal_polynomial, candidate_minimal_polynomial))
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
    struct Number*a, struct RationalPolynomial*a_minimal_polynomial)
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
        struct RationalPolynomial*radicand_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a->left);
        struct Number**radicand_conjugates =
            number_conjugates(local_stack, output_stack, a->left, radicand_minimal_polynomial);
        for (size_t i = 0; i < radicand_minimal_polynomial->coefficient_count - 1; ++i)
        {
            radicand_conjugates[i] = number_exponentiate(local_stack, output_stack,
                radicand_conjugates[i], &a->right->value);
        }
        size_t roots_of_unity_count = integer_to_size_t(a->right->value.denominator);
        struct Number**roots_of_unity =
            get_roots_of_unity(output_stack, local_stack, a->right->value.denominator);
        struct Number**out = number_sum_or_product_conjugates(number_multiply, output_stack,
            local_stack, radicand_conjugates, radicand_minimal_polynomial->coefficient_count - 1,
            roots_of_unity, roots_of_unity_count, a_minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalPolynomial*left_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a->left);
        struct Number**left_conjugates =
            number_conjugates(local_stack, output_stack, a->left, left_minimal_polynomial);
        struct RationalPolynomial*right_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a->right);
        struct Number**right_conjugates =
            number_conjugates(local_stack, output_stack, a->right, right_minimal_polynomial);
        struct Number**out = number_sum_or_product_conjugates(number_multiply, output_stack,
            local_stack, left_conjugates, left_minimal_polynomial->coefficient_count - 1,
            right_conjugates, right_minimal_polynomial->coefficient_count - 1,
            a_minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '+':
    {
        if (a->term_count == 1)
        {
            return number_conjugates(output_stack, local_stack, a->terms[0], a_minimal_polynomial);
        }
        void*local_stack_savepoint = local_stack->cursor;
        --a->term_count;
        struct RationalPolynomial*left_terms_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a);
        struct Number**left_terms_conjugates =
            number_conjugates(local_stack, output_stack, a, left_terms_minimal_polynomial);
        struct RationalPolynomial*right_term_minimal_polynomial =
            number_minimal_polynomial(local_stack, output_stack, a->terms[a->term_count]);
        struct Number**right_term_conjugates = number_conjugates(local_stack, output_stack,
            a->terms[a->term_count], right_term_minimal_polynomial);
        ++a->term_count;
        struct Number**out = number_sum_or_product_conjugates(number_add, output_stack, local_stack,
            left_terms_conjugates, left_terms_minimal_polynomial->coefficient_count - 1,
            right_term_conjugates, right_term_minimal_polynomial->coefficient_count - 1,
            a_minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}