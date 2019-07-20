#include "declarations.h"

void number_calculate_minimal_polynomial_from_annulling_polynomial(struct PoolSet*pool_set,
    struct Stack*stack_a, struct Stack*stack_b, struct RationalPolynomial*annulling_polynomial,
    struct Number*a)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct RationalPolynomial**candidates = stack_slot_allocate(stack_a,
        annulling_polynomial->coefficient_count * sizeof(struct RationalPolynomial*),
        _Alignof(struct RationalPolynomial*));
    size_t candidate_count =
        rational_polynomial_factor(stack_a, stack_b, annulling_polynomial, candidates);
    if (candidate_count == 1)
    {
        a->minimal_polynomial = rational_polynomial_copy_to_number_memory(pool_set, candidates[0]);
        stack_a->cursor = stack_a_savepoint;
        return;
    }
    struct FloatInterval zero_estimate = { &float_zero, &float_zero };
    struct RectangularEstimate zero_rectangular_estimate = { &zero_estimate, &zero_estimate };
    struct Rational*interval_size = &rational_one;
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct RectangularEstimate evaluation_estimate;
            rational_polynomial_estimate_evaluation(pool_set, stack_a, stack_b,
                &evaluation_estimate, candidates[i], a, interval_size);
            if (rectangular_estimates_are_disjoint(stack_a, stack_b, &evaluation_estimate,
                &zero_rectangular_estimate))
            {
                candidate_count -= 1;
                candidates[i] = candidates[candidate_count];
                if (candidate_count == 1)
                {
                    a->minimal_polynomial =
                        rational_polynomial_copy_to_number_memory(pool_set, candidates[0]);
                    stack_a->cursor = stack_a_savepoint;
                    return;
                }
            }
            else
            {
                ++i;
            }
        }
        interval_size->denominator = integer_doubled(stack_a, interval_size->denominator);
    }
}

void number_calculate_minimal_polynomial_from_algebraic_form(struct PoolSet*pool_set,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a,
    struct AlgebraicNumber*algebraic_form,
    struct RationalPolynomial*generator_minimal_polynomials[2])
{
    void*stack_a_savepoint = stack_a->cursor;
    size_t term_count_upper_bound = (generator_minimal_polynomials[0]->coefficient_count - 1) *
        (generator_minimal_polynomials[1]->coefficient_count - 1) + 1;
    struct Matrix matrix;
    matrix.rows = stack_slot_allocate(stack_a, term_count_upper_bound * sizeof(struct Rational**),
        _Alignof(struct Rational**));
    matrix.height = 0;
    size_t(**terms_present)[2] = stack_slot_allocate(stack_a,
        term_count_upper_bound * sizeof(size_t(*)[2]), _Alignof(size_t(*)[2]));
    size_t terms_present_count = 1;
    struct AlgebraicNumber*power = STACK_SLOT_ALLOCATE(stack_a, struct AlgebraicNumber);
    power->term_coefficient = &rational_one;
    power->generator_degrees[0] = 0;
    power->generator_degrees[1] = 0;
    power->next_term = 0;
    struct AlgebraicNumber**powers =
        stack_slot_allocate(stack_a, term_count_upper_bound * sizeof(struct AlgebraicNumber*),
            _Alignof(struct AlgebraicNumber*));
    bool constant_is_present = false;
    while (!constant_is_present || matrix.height < terms_present_count)
    {
        power = algebraic_number_multiply(stack_a, stack_b, power, algebraic_form,
            generator_minimal_polynomials);
        powers[matrix.height] = power;
        struct AlgebraicNumber*power_term = powers[matrix.height];
        matrix.rows[matrix.height] = stack_slot_allocate(stack_a,
            term_count_upper_bound * sizeof(struct Rational*), _Alignof(struct Rational*));
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
    struct RationalPolynomial**augmentation = stack_slot_allocate(stack_a,
        matrix.height * sizeof(struct RationalPolynomial*), _Alignof(struct RationalPolynomial*));
    for (size_t i = 0; i < matrix.height; ++i)
    {
        array_reverse(matrix.rows[i], matrix.width);
        augmentation[i] = polynomial_allocate(stack_a, i + 2);
        for (size_t j = 0; j <= i; ++j)
        {
            augmentation[i]->coefficients[j] = &rational_zero;
        }
        augmentation[i]->coefficients[i + 1] = &rational_one;
    }
    matrix_row_echelon_form(stack_a, stack_b, &matrix, augmentation);
    struct RationalPolynomial*annulling_polynomial = augmentation[matrix.height - 1];
    annulling_polynomial->coefficients[0] = rational_subtract(stack_a, stack_b,
        annulling_polynomial->coefficients[0], matrix.rows[matrix.height - 1][matrix.width - 1]);
    number_calculate_minimal_polynomial_from_annulling_polynomial(pool_set, stack_b, stack_a,
        annulling_polynomial, a);
    stack_a->cursor = stack_a_savepoint;
}

void number_calculate_minimal_polynomial(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a)
{
    if (a->minimal_polynomial)
    {
        return;
    }
    switch (a->operation)
    {
    case 'r':
        a->minimal_polynomial = polynomial_allocate(&polynomial_stack, 2);
        a->minimal_polynomial->coefficients[0]->numerator =
            integer_copy_to_pool(pool_set, a->value.numerator);
        a->minimal_polynomial->coefficients[0]->numerator->sign =
            -a->minimal_polynomial->coefficients[0]->numerator->sign;
        a->minimal_polynomial->coefficients[0]->denominator =
            integer_copy_to_pool(pool_set, a->value.denominator);
        a->minimal_polynomial->coefficients[1]->numerator = pool_integer_initialize(pool_set, 1, 1);
        a->minimal_polynomial->coefficients[1]->denominator =
            a->minimal_polynomial->coefficients[1]->numerator;
        return;
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a->left);
        size_t surd_index = integer_to_size_t(a->left->value.denominator);
        struct RationalPolynomial*annulling_polynomial = polynomial_allocate(stack_a,
            surd_index * (a->left->minimal_polynomial->coefficient_count - 1) + 1);
        annulling_polynomial->coefficients[0] = a->left->minimal_polynomial->coefficients[0];
        for (size_t i = 1; i < a->left->minimal_polynomial->coefficient_count; ++i)
        {
            size_t annulling_polynomial_coefficient_index = i * surd_index;
            annulling_polynomial->coefficients[annulling_polynomial_coefficient_index] =
                a->left->minimal_polynomial->coefficients[i];
            for (size_t j = 1; j < surd_index; ++j)
            {
                annulling_polynomial->coefficients[annulling_polynomial_coefficient_index - j] =
                    &rational_zero;
            }
        }
        number_calculate_minimal_polynomial_from_annulling_polynomial(pool_set, stack_a,
            stack_b, annulling_polynomial, a);
        stack_a->cursor = stack_a_savepoint;
        return;
    }
    case '*':
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct AlgebraicNumber*algebraic_form =
            STACK_SLOT_ALLOCATE(stack_a, struct AlgebraicNumber);
        algebraic_form->term_coefficient = &rational_one;
        algebraic_form->generator_degrees[0] = 1;
        algebraic_form->generator_degrees[1] = 1;
        algebraic_form->next_term = 0;
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a->left);
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a->right);
        number_calculate_minimal_polynomial_from_algebraic_form(pool_set, stack_a, stack_b, a,
            algebraic_form, (struct RationalPolynomial*[2]) { a->left->minimal_polynomial,
                a->right->minimal_polynomial });
        stack_a->cursor = stack_a_savepoint;
        return;
    }
    case '+':
        crash("number_minimal_polynomial case not yet implemented.");
    }
}

void number_calculate_sum_or_product_conjugates(struct Number*(operation)(struct PoolSet*,
    struct Stack*, struct Stack*, struct Number*, struct Number*),
    struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b, struct Number*a)
{
    number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a);
    number_calculate_conjugates(pool_set, stack_a, stack_b, a->left);
    number_calculate_conjugates(pool_set, stack_a, stack_b, a->right);
    struct Number*conjugate = a;
    struct Number*left_conjugate = a->left;
    struct Number*right_conjugate = a->right->next;
    size_t conjugate_count = 1;
    while (left_conjugate)
    {
        while (right_conjugate)
        {
            struct Number*candidate = operation(pool_set, stack_a, stack_b,
                number_copy(pool_set, left_conjugate), number_copy(pool_set, right_conjugate));
            number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, candidate);
            if (!rational_polynomial_equals(a->minimal_polynomial, candidate->minimal_polynomial))
            {
                conjugate->next = candidate;
                ++conjugate_count;
                if (conjugate_count == a->minimal_polynomial->coefficient_count - 1)
                {
                    return;
                }
                conjugate = conjugate->next;
            }
            else
            {
                number_free(pool_set, candidate);
            }
            right_conjugate = right_conjugate->next;
        }
        left_conjugate = left_conjugate->next;
        right_conjugate = a->right;
    }
    crash("Not enough conjugates found.");
}

void number_calculate_conjugates(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a)
{
    if (a->next)
    {
        return;
    }
    switch (a->operation)
    {
    case 'r':
        return;
    case '^':
    {
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a);
        number_calculate_conjugates(pool_set, stack_a, stack_b, a->left);
        size_t roots_of_unity_count = integer_to_size_t(a->right->value.denominator);
        struct Number**roots_of_unity =
            get_roots_of_unity(stack_a, stack_b, a->right->value.denominator);
        struct Number*radicand_conjugate = a->left->next;
        struct Number*conjugate = a;
        size_t conjugate_count = 1;
        for (size_t i = 0; i < roots_of_unity_count; ++i)
        {
            while (radicand_conjugate)
            {
                struct Number*surd = number_allocate(pool_set);
                surd->operation = '^';
                surd->left = number_copy(pool_set, radicand_conjugate);
                surd->right = number_copy(pool_set, a->left);
                struct Number*candidate = number_multiply(pool_set, stack_a, stack_b, surd,
                    number_copy(pool_set, radicand_conjugate));
                number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, candidate);
                if (!rational_polynomial_equals(a->minimal_polynomial,
                    candidate->minimal_polynomial))
                {
                    conjugate->next = candidate;
                    ++conjugate_count;
                    if (conjugate_count == a->minimal_polynomial->coefficient_count - 1)
                    {
                        return;
                    }
                    conjugate = conjugate->next;
                }
                else
                {
                    number_free(pool_set, candidate);
                }
                radicand_conjugate = radicand_conjugate->next;
            }
            radicand_conjugate = a->left;
        }
        crash("Not enough conjugates found.");
    }
    case '*':
    {
        number_calculate_sum_or_product_conjugates(number_multiply, pool_set, stack_a, stack_b, a);
        return;
    }
    case '+':
    {
        number_calculate_sum_or_product_conjugates(number_add, pool_set, stack_a, stack_b, a);
    }
    }
}