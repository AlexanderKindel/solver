#include "declarations.h"

struct Number*number_copy(struct Stack*output_stack, struct Number*a)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = a->operation;
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

struct RationalPolynomial*number_a_in_terms_of_b(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*a_minimal_polynomial,
    struct Number*b, struct RationalPolynomial*b_minimal_polynomial)
{
    if (a->operation == 'r')
    {
        struct RationalPolynomial*out = polynomial_allocate(output_stack, 1);
        out->coefficients[0] = rational_copy(output_stack, &a->value);
        return out;
    }
    if ((b_minimal_polynomial->coefficient_count - 1) %
        (a_minimal_polynomial->coefficient_count - 1) != 0)
    {
        return 0;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*nested_minimal_polynomial =
        rational_polynomial_to_nested_polynomial(local_stack, a_minimal_polynomial);
    struct NestedPolynomial**minimal_polynomial_factors = stack_slot_allocate(local_stack,
        (nested_minimal_polynomial->coefficient_count - 1) * sizeof(struct NestedPolynomial*),
        _Alignof(struct NestedPolynomial*));
    size_t candidate_factor_count = number_field_polynomial_factor(local_stack, output_stack,
        nested_minimal_polynomial, b_minimal_polynomial, minimal_polynomial_factors);
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
    struct Number**conjugates =
        number_conjugates(output_stack, local_stack, a, a_minimal_polynomial);
    size_t conjugate_count = a_minimal_polynomial->coefficient_count - 1;
    size_t*conjugate_indices = ARRAY_ALLOCATE(local_stack, conjugate_count, size_t*);
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

struct Number*number_evaluate(struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    if (a->operation == 'r')
    {
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