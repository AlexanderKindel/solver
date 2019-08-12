#include "declarations.h"

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
                    rational_integer_divide(local_stack, output_stack, pi.min, product_index);
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
                &(struct Rational){ &one, product_index });
            struct RationalInterval product_radicand_rational_argument_estimate;
            number_rational_argument_estimate(local_stack, output_stack,
                &product_radicand_rational_argument_estimate, product_radicand, pi.min);
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