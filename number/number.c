#include "declarations.h"

struct Number*number_copy(struct Stack*output_stack, struct Number*a)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = a->operation;
    if (a->minimal_polynomial)
    {
        out->minimal_polynomial = rational_polynomial_copy(output_stack, a->minimal_polynomial);
    }
    else
    {
        out->minimal_polynomial = 0;
    }
    switch (a->operation)
    {
    case 'r':
        out->value = rational_copy(output_stack, &a->value);
        break;
    case '^':
        out->radicand = number_copy(output_stack, a->radicand);
        out->index = integer_copy(output_stack, a->index);
        break;
    case '+':
        if (a->generator)
        {
            out->generator = number_copy(output_stack, a->generator);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, a->element_count, struct RationalPolynomial*);
            for (size_t i = 0; i < a->element_count; ++i)
            {
                out->terms_in_terms_of_generator[i] =
                    rational_polynomial_copy(output_stack, a->terms_in_terms_of_generator[i]);
            }
        }
        else
        {
            out->generator = 0;
        }
    case '*':
        out->element_count = a->element_count;
        out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
        for (size_t i = 0; i < a->element_count; ++i)
        {
            out->elements[i] = number_copy(output_stack, a->elements[i]);
        }
    }
    return out;
}

struct Number*number_rational_initialize(struct Stack*output_stack, struct Rational*value)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = 'r';
    out->value = rational_copy(output_stack, value);
    out->minimal_polynomial = POLYNOMIAL_ALLOCATE(output_stack, 2, struct Rational);
    out->minimal_polynomial->coefficients[0].numerator =
        integer_negate(output_stack, value->numerator);
    out->minimal_polynomial->coefficients[0].denominator =
        integer_copy(output_stack, value->denominator);
    out->minimal_polynomial->coefficients[1] = rational_one;
    return out;
}

int8_t sum_or_product_formal_compare(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct Number*a, struct Number*b)
{
    if (a->element_count > b->element_count)
    {
        return 1;
    }
    if (a->element_count < b->element_count)
    {
        return -1;
    }
    for (size_t i = 0; i < a->element_count; ++i)
    {
        int8_t out =
            number_formal_compare(local_stack_a, local_stack_b, a->elements[i], b->elements[i]);
        if (out)
        {
            return out;
        }
    }
    return 0;
}

int8_t number_formal_compare(struct Stack*restrict local_stack_a,
    struct Stack*restrict local_stack_b, struct Number*a, struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (b->operation == 'r')
        {
            return rational_compare(local_stack_a, local_stack_b, &a->value, &b->value);
        }
        return -1;
    case '^':
        switch (b->operation)
        {
        case 'r':
            return 1;
        case '^':
        {
            int8_t out = integer_compare(local_stack_a, local_stack_b, a->index, b->index);
            if (out)
            {
                return -out;
            }
            return number_formal_compare(local_stack_a, local_stack_b, a->radicand, b->radicand);
        }
        }
        return -1;
    case '*':
        switch (b->operation)
        {
        case '+':
            return -1;
        case '*':
            return sum_or_product_formal_compare(local_stack_a, local_stack_b, a, b);
        }
        return 1;
    case '+':
        if (b->operation != '+')
        {
            return 1;
        }
        return sum_or_product_formal_compare(local_stack_a, local_stack_b, a, b);
    }
    crash("Number operation not recognized.");
}

struct RationalPolynomial*number_annulling_polynomial_to_minimal_polynomial(
    struct Stack*restrict output_stack, struct Stack*restrict local_stack,
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
    struct FloatInterval zero_estimate = { float_zero, float_zero };
    struct Region a_estimate =
        number_get_region_estimate(local_stack, output_stack, a, &rational_one);
    struct Rational interval_size_for_evaluation =
    { &one, size_t_to_integer(local_stack, annulling_polynomial->coefficient_count) };
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct Region evaluation_estimate =
                rational_polynomial_evaluate_at_region_estimate(local_stack, output_stack,
                    candidates[i], &a_estimate, &interval_size_for_evaluation);
            if (regions_are_disjoint(output_stack, local_stack, &evaluation_estimate,
                &(struct Region) { zero_estimate, zero_estimate }))
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
        number_halve_region_estimate_dimensions(local_stack, output_stack, &a_estimate, a);
        interval_size_for_evaluation.denominator =
            integer_double(local_stack, interval_size_for_evaluation.denominator);
    }
}

size_t trailing_factor_to_string(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number**factors, size_t factor_count)
{
    size_t char_count = 0;
    for (size_t i = 0; i < factor_count; ++i)
    {
        *(char*)ALLOCATE(output_stack, char) = '*';
        char_count += 1 + number_to_string(output_stack, local_stack, factors[i]);
    }
    return char_count;
}

size_t product_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number**factors, size_t factor_count)
{
    size_t char_count = number_to_string(output_stack, local_stack, factors[0]);
    return char_count +
        trailing_factor_to_string(output_stack, local_stack, factors + 1, factor_count - 1);
}

size_t number_to_string(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
    {
        size_t char_count = integer_to_string(output_stack, local_stack, a->value.numerator);
        if (!integer_equals(a->value.denominator, &one))
        {
            *(char*)ALLOCATE(output_stack, char) = '/';
            char_count += 1 + integer_to_string(output_stack, local_stack, a->value.denominator);
        }
        return char_count;
    }
    case '^':
    {
        size_t char_count;
        if (a->radicand->operation != 'r' || a->radicand->value.numerator->sign < 0)
        {
            *(char*)ALLOCATE(output_stack, char) = '(';
            char_count = 2 + number_to_string(output_stack, local_stack, a->radicand);
            *(char*)ALLOCATE(output_stack, char) = ')';
        }
        else
        {
            char_count = number_to_string(output_stack, local_stack, a->radicand);
        }
        char*exponent = ARRAY_ALLOCATE(output_stack, 4, char);
        exponent[0] = '^';
        exponent[1] = '(';
        exponent[2] = '1';
        exponent[3] = '/';
        char_count += 5 + integer_to_string(output_stack, local_stack, a->index);
        *(char*)ALLOCATE(output_stack, char) = ')';
        return char_count;
    }
    case '*':
    {
        size_t char_count;
        if (a->elements[0]->operation == 'r')
        {
            if (integer_equals(a->elements[0]->value.numerator, INT(1, -1)))
            {
                *(char*)ALLOCATE(output_stack, char) = '-';
                char_count = 1 + product_to_string(output_stack, local_stack, a->elements + 1,
                    a->element_count - 1);
            }
            else if (integer_equals(a->elements[0]->value.numerator, &one))
            {
                char_count = product_to_string(output_stack, local_stack, a->elements + 1,
                    a->element_count - 1);
            }
            else
            {
                char_count = integer_to_string(output_stack, local_stack,
                    a->elements[0]->value.numerator);
                char*next_factor = local_stack->cursor;
                size_t next_factor_char_count =
                    number_to_string(local_stack, output_stack, a->elements[1]);
                if (next_factor[0] != '(')
                {
                    *(char*)ALLOCATE(output_stack, char) = '*';
                    char_count += 1;
                }
                char*next_factor_out = ARRAY_ALLOCATE(output_stack, next_factor_char_count, char);
                memcpy(next_factor_out, next_factor, next_factor_char_count);
                char_count += next_factor_char_count + trailing_factor_to_string(output_stack,
                    local_stack, a->elements + 2, a->element_count - 2);
                local_stack->cursor = next_factor;
            }
            if (!integer_equals(a->elements[0]->value.denominator, &one))
            {
                *(char*)ALLOCATE(output_stack, char) = '/';
                char_count += 1 + integer_to_string(output_stack, local_stack,
                    a->elements[0]->value.denominator);
            }
        }
        else
        {
            char_count =
                product_to_string(output_stack, local_stack, a->elements, a->element_count);
        }
        return char_count;
    }
    case '+':
    {
        size_t char_count = number_to_string(output_stack, local_stack, a->elements[0]);
        for (size_t i = 1; i < a->element_count; ++i)
        {
            char*sum_string = local_stack->cursor;
            size_t sum_char_count = number_to_string(local_stack, output_stack, a->elements[i]);
            if (sum_string[0] != '-')
            {
                *(char*)ALLOCATE(output_stack, char) = '+';
                char_count += 1;
            }
            char*sum_string_out = ARRAY_ALLOCATE(output_stack, sum_char_count, char);
            memcpy(sum_string_out, sum_string, sum_char_count);
            char_count += sum_char_count;
            local_stack->cursor = sum_string;
        }
        return char_count;
    }
    default:
        crash("Number operation not recognized.");
    }
}

#include "number/number_add/number_add.c"
#include "number/number_divide/number_divide.c"
#include "number/number_estimate/number_estimate.c"
#include "number/number_exponentiate/number_exponentiate.c"
#include "number/number_multiply/number_multiply.c"
#include "number/roots_of_unity/roots_of_unity.c"