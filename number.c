#include "declarations.h"

struct Number*number_allocate(struct PoolSet*pool_set)
{
    struct Number*out = pool_slot_allocate(pool_set, sizeof(struct Number));
    memset(out, 0, sizeof(struct Number));
    return out;
}

void number_node_free_during_parse(struct PoolSet*pool_set, struct Number*a)
{
    if (a->operation == 'r')
    {
        rational_free(pool_set, &a->value);
    }
    pool_slot_free(pool_set, a, sizeof(struct Number));
}

void number_node_free(struct PoolSet*pool_set, struct Number*a)
{
    if (a->operation != 'r')
    {
        float_interval_free(pool_set, a->real_part_estimate);
        float_interval_free(pool_set, a->imaginary_part_estimate);
    }
    float_interval_free(pool_set, a->argument_estimate);
    float_interval_free(pool_set, a->magnitude_estimate);
    if (a->next)
    {
        number_free(pool_set, a->next);
    }
    number_node_free_during_parse(pool_set, a);
}

void number_free(struct PoolSet*pool_set, struct Number*a)
{
    if (a->operation != 'r')
    {
        number_free(pool_set, a->left);
        number_free(pool_set, a->right);
    }
    number_node_free(pool_set, a);
}

struct Number*number_copy(struct PoolSet*pool_set, struct Number*a)
{
    struct Number*out = number_allocate(pool_set);
    memcpy(out, a, sizeof(struct Number));
    if (a->operation == 'r')
    {
        rational_move_to_pool(pool_set, &out->value);
    }
    else
    {
        out->left = number_copy(pool_set, a->left);
        out->right = number_copy(pool_set, a->right);
        if (a->real_part_estimate)
        {
            float_interval_move_to_pool(pool_set, out->real_part_estimate);
        }
        if (a->imaginary_part_estimate)
        {
            float_interval_move_to_pool(pool_set, out->imaginary_part_estimate);
        }
    }
    if (a->argument_estimate)
    {
        float_interval_move_to_pool(pool_set, out->argument_estimate);
    }
    if (a->magnitude_estimate)
    {
        float_interval_move_to_pool(pool_set, out->magnitude_estimate);
    }
    return out;
}

struct RationalPolynomial*number_a_in_terms_of_b(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct Number*b)
{
    if (a->operation == 'r')
    {
        struct RationalPolynomial*out = polynomial_allocate(output_stack, 1);
        out->coefficients[0] = rational_copy_to_stack(output_stack, &a->value);
        return out;
    }
    number_calculate_minimal_polynomial(pool_set, output_stack, local_stack, a);
    number_calculate_minimal_polynomial(pool_set, output_stack, local_stack, b);
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
    number_calculate_conjugates(pool_set, output_stack, local_stack, a);
    size_t conjugate_count = 0;
    struct Number**conjugates = STACK_SLOT_ALLOCATE(local_stack, struct Number*);
    struct Number*conjugate = a->next;
    while (conjugate)
    {
        conjugates[conjugate_count] = conjugate;
        STACK_SLOT_ALLOCATE(local_stack, struct Number*);
        ++conjugate_count;
        conjugate = conjugate->next;
    }
    size_t*conjugate_indices =
        stack_slot_allocate(local_stack, conjugate_count * sizeof(size_t*), _Alignof(size_t*));
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
            number_rectangular_estimate(pool_set, local_stack, output_stack, &a_estimate, a,
                interval_size);
            struct RectangularEstimate factor_at_b_estimate;
            rational_polynomial_estimate_evaluation(pool_set, local_stack, output_stack,
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
                number_rectangular_estimate(pool_set, local_stack, output_stack,
                    &candidate_estimate, conjugate, interval_size);
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

size_t number_string(struct Stack*output_stack, struct Stack*local_stack, struct Number*number)
{
    switch (number->operation)
    {
    case 'r':
    {
        size_t char_count = integer_string(output_stack, local_stack, number->value.numerator);
        if (!integer_equals(number->value.denominator, &one))
        {
            *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = '/';
            char_count += 1 + integer_string(output_stack, local_stack, number->value.denominator);
        }
        return char_count;
    }
    case '^':
    {
        size_t char_count;
        if (number->left->operation != 'r' || number->left->value.numerator->sign < 0)
        {
            *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = '(';
            char_count = 2 + number_string(output_stack, local_stack, number->left);
            *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = ')';
        }
        else
        {
            char_count = number_string(output_stack, local_stack, number->left);
        }
        char*exponent = stack_slot_allocate(output_stack, 2 * sizeof(char), _Alignof(char));
        exponent[0] = '^';
        exponent[1] = '(';
        char_count += 3 + number_string(output_stack, local_stack, number->right);
        *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = ')';
        return char_count;
    }
    case '*':
    {
        size_t char_count;
        if (number->left->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            if (integer_equals(number->left->value.numerator,
                stack_integer_initialize(local_stack, 1, -1)))
            {
                *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = '-';
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
                    STACK_SLOT_ALLOCATE(output_stack, char);
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
                *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = '/';
                char_count +=
                    1 + integer_string(output_stack, local_stack, number->left->value.denominator);
            }
            local_stack->cursor = local_stack_savepoint;
        }
        else
        {
            char_count = number_string(output_stack, local_stack, number->left) + 1;
            *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = '*';
            char_count += number_string(output_stack, local_stack, number->right);
        }
        return char_count;
    }
    case '+':
    {
        size_t char_count = number_string(output_stack, local_stack, number->left) + 1;
        *(char*)STACK_SLOT_ALLOCATE(output_stack, char) = '+';
        return char_count + number_string(output_stack, local_stack, number->right);
    }
    default:
        crash("Number operation not recognized.");
    }
}