#include "declarations.h"

struct RationalPolynomial*product_minimal_polynomial(struct Stack*output_stack,
    struct Stack*local_stack, struct Number*a, struct RationalPolynomial*left_minimal_polynomial,
    struct RationalPolynomial*right_minimal_polynomial)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct NestedPolynomial*t =
        polynomial_allocate(local_stack, left_minimal_polynomial->coefficient_count);
    for (size_t i = 0; i < left_minimal_polynomial->coefficient_count; ++i)
    {
        t->coefficients[i] =
            polynomial_allocate(local_stack, left_minimal_polynomial->coefficient_count - i);
        size_t degree = t->coefficients[i]->coefficient_count - 1;
        for (size_t j = 0; j < degree; ++j)
        {
            t->coefficients[i]->coefficients[j] = &rational_zero;
        }
        t->coefficients[i]->coefficients[degree] =
            left_minimal_polynomial->coefficients[degree - i];
    }
    struct RationalPolynomial*out =
        number_minimal_polynomial_from_annulling_polynomial(output_stack, local_stack,
            nested_polynomial_resultant(local_stack, output_stack, t,
                rational_polynomial_to_nested_polynomial(local_stack, right_minimal_polynomial)),
            a);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_product_initialize(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*left_factor, struct Number*right_factor)
{
    struct Number*out = ALLOCATE(output_stack, struct Number);
    out->operation = '*';
    out->element_count = 2;
    out->elements = ARRAY_ALLOCATE(output_stack, 2, struct Number*);
    out->elements[0] = number_copy(output_stack, left_factor);
    out->elements[1] = number_copy(output_stack, right_factor);
    if (number_formal_compare(output_stack, local_stack, out->elements[0], out->elements[1]) > 0)
    {
        POINTER_SWAP(out->elements[0], out->elements[1]);
    }
    out->minimal_polynomial = product_minimal_polynomial(output_stack, local_stack, out,
        left_factor->minimal_polynomial, right_factor->minimal_polynomial);
    return out;
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
            rational_multiply(local_stack, output_stack, a->value, b));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '^':
        out = ALLOCATE(output_stack, struct Number);
        out->operation = '*';
        out->element_count = 2;
        out->elements = ARRAY_ALLOCATE(output_stack, 2, struct Number*);
        out->elements[0] = number_rational_initialize(output_stack, b);
        out->elements[1] = number_copy(output_stack, a);
        break;
    case '*':
        size_t factors_to_copy_count;
        if (a->elements[0]->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Rational*coefficient =
                rational_multiply(local_stack, output_stack, a->elements[0]->value, b);
            if (rational_equals(coefficient, &rational_one))
            {
                if (a->element_count == 2)
                {
                    out = number_copy(output_stack, a->elements[1]);
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
                out = ALLOCATE(output_stack, struct Number);
                out->element_count = a->element_count - 1;
                out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
            }
            else
            {
                out = ALLOCATE(output_stack, struct Number);
                out->element_count = a->element_count;
                out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
                out->elements[0] = number_rational_initialize(output_stack, coefficient);
            }
            factors_to_copy_count = a->element_count - 1;
            local_stack->cursor = local_stack_savepoint;
        }
        else
        {
            out = ALLOCATE(output_stack, struct Number);
            out->element_count = a->element_count + 1;
            out->elements = ARRAY_ALLOCATE(output_stack, out->element_count, struct Number*);
            out->elements[0] = number_rational_initialize(output_stack, b);
            factors_to_copy_count = a->element_count;
        }
        out->operation = '*';
        for (size_t i = 1; i <= factors_to_copy_count; ++i)
        {
            out->elements[out->element_count - i] =
                number_copy(output_stack, a->elements[a->element_count - i]);
        }
        break;
    default:
        out = ALLOCATE(output_stack, struct Number);
        out->operation = a->operation;
        if (a->operation == '+')
        {
            out->generator = number_copy(output_stack, a->generator);
            out->terms_in_terms_of_generator =
                ARRAY_ALLOCATE(output_stack, a->element_count, struct RationalPolynomial*);
            for (size_t i = 0; i < a->element_count; ++i)
            {
                out->terms_in_terms_of_generator[i] =
                    rational_polynomial_rational_multiply(output_stack, local_stack,
                        a->terms_in_terms_of_generator[i], b);
            }
        }
        out->element_count = a->element_count;
        out->elements = ARRAY_ALLOCATE(output_stack, a->element_count, struct Number*);
        for (size_t i = 0; i < a->element_count; ++i)
        {
            out->elements[i] =
                number_rational_multiply(output_stack, local_stack, a->elements[i], b);
        }
    }
    void*local_stack_savepoint = local_stack->cursor;
    out->minimal_polynomial =
        polynomial_allocate(local_stack, a->minimal_polynomial->coefficient_count);
    struct Rational*power = &rational_one;
    for (size_t i = 0; i < a->minimal_polynomial->coefficient_count; ++i)
    {
        out->minimal_polynomial->coefficients[i] = rational_divide(local_stack, output_stack,
            a->minimal_polynomial->coefficients[i], power);
        power = rational_unreduced_multiply(local_stack, output_stack, power, b, 0);
    }
    out->minimal_polynomial = rational_polynomial_rational_multiply(output_stack, local_stack,
        out->minimal_polynomial, rational_reciprocal(local_stack,
            out->minimal_polynomial->coefficients[out->minimal_polynomial->coefficient_count - 1]));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*factor_consolidate(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a, struct Number*b)
{
    if (a->operation == 'r')
    {
        if (b->operation == 'r')
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*out = number_rational_initialize(output_stack,
                rational_multiply(local_stack, output_stack, a->value, b->value));
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        return 0;
    }
    if (b->operation == 'r')
    {
        return 0;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct ExtendedGCDInfo info;
    integer_extended_gcd(local_stack, output_stack, &info, a->index, b->index);
    struct Integer*product_index = integer_euclidean_quotient(local_stack, output_stack,
        integer_multiply(local_stack, output_stack, a->index, b->index), info.gcd);
    struct Number*radicand_factor_a = number_integer_exponentiate(local_stack, output_stack,
        a->radicand, integer_magnitude(local_stack, info.b_over_gcd));
    struct Number*radicand_factor_b = number_integer_exponentiate(local_stack, output_stack,
        b->radicand, integer_magnitude(local_stack, info.a_over_gcd));
    struct RationalInterval*radicand_factor_a_argument_estimate =
        number_rational_argument_estimate(local_stack, output_stack, radicand_factor_a,
            &rational_one);
    struct RationalInterval*radicand_factor_b_argument_estimate =
        number_rational_argument_estimate(local_stack, output_stack, radicand_factor_b,
            &rational_one);
    if (!integer_equals(a->index, b->index))
    {
        struct Rational*radicand_factor_argument_interval_size =
            rational_integer_divide(local_stack, output_stack, pi.min, product_index);
        if (rational_compare(output_stack, local_stack,
            rational_integer_divide(local_stack, output_stack,
                radicand_factor_a_argument_estimate->max, product_index),
            number_rational_argument_estimate(local_stack, output_stack, a,
                radicand_factor_argument_interval_size)->min) < 0 ||
            rational_compare(output_stack, local_stack,
                rational_integer_divide(local_stack, output_stack,
                    radicand_factor_b_argument_estimate->max, product_index),
                number_rational_argument_estimate(local_stack, output_stack, b,
                    radicand_factor_argument_interval_size)->min) < 0)
        {
            local_stack->cursor = local_stack_savepoint;
            return 0;
        }
    }
    struct Number*product_radicand =
        number_multiply(local_stack, output_stack, radicand_factor_a, radicand_factor_b);
    struct Number*out = number_root(local_stack, output_stack, product_radicand, product_index);
    if (rational_compare(output_stack, local_stack,
        number_rational_argument_estimate(local_stack, output_stack, product_radicand, pi.min)->max,
        rational_add(local_stack, output_stack, radicand_factor_a_argument_estimate->min,
            radicand_factor_b_argument_estimate->min)) < 0)
    {
        out = number_multiply(output_stack, local_stack, out,
            get_roots_of_unity(output_stack, local_stack, product_index)[1]);
    }
    else
    {
        out = number_copy(output_stack, out);
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_multiply(struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Number*b)
{
    if (b->operation == 'r')
    {
        return number_rational_multiply(output_stack, local_stack, a, b->value);
    }
    switch (a->operation)
    {
    case 'r':
    {
        return number_rational_multiply(output_stack, local_stack, b, a->value);
    }
    case '^':
    {
        if (b->operation == '^')
        {
            void*output_stack_savepoint = output_stack->cursor;
            struct Number*out = factor_consolidate(output_stack, local_stack, a, b);
            if (out)
            {
                return out;
            }
            output_stack->cursor = output_stack_savepoint;
            if (number_formal_compare(output_stack, local_stack, a, b) < 0)
            {
                return number_product_initialize(output_stack, local_stack, a, b);
            }
            return number_product_initialize(output_stack, local_stack, b, a);
        }
        break;
    }
    case '*':
        switch (b->operation)
        {
        case '^':
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number**factors = ARRAY_ALLOCATE(local_stack, a->element_count, struct Number*);
            memcpy(factors, a->elements, a->element_count * sizeof(struct Number*));
            size_t factor_count = a->element_count + 1;
            size_t factor_index = 0;
            struct Number*consolidated_factors = b;
            while (factor_index < a->element_count)
            {
                if (a->elements[factor_index])
                {
                    struct Number*consolidation_attempt =
                        factor_consolidate(local_stack, output_stack, a->elements[factor_index], b);
                    if (consolidation_attempt)
                    {
                        consolidated_factors = consolidation_attempt;
                        a->elements[factor_index] = 0;
                        --factor_count;
                        factor_index = 0;
                        continue;
                    }
                }
                ++factor_index;
            }
            if (factor_count == 1)
            {
                consolidated_factors = number_copy(output_stack, consolidated_factors);
                local_stack->cursor = local_stack_savepoint;
                return consolidated_factors;
            }
            struct Number*out = ALLOCATE(output_stack, struct Number);
            out->operation = '*';
            out->elements = ARRAY_ALLOCATE(output_stack, factor_count, struct Number*);
            out->element_count = 0;
            factor_index = 0;
            while (true)
            {
                if (factors[factor_index])
                {
                    if (number_formal_compare(output_stack, local_stack, factors[factor_index],
                        consolidated_factors) < 0)
                    {
                        out->elements[out->element_count] =
                            number_copy(output_stack, factors[factor_index]);
                        ++out->element_count;
                    }
                    else
                    {
                        out->elements[out->element_count] =
                            number_copy(output_stack, consolidated_factors);
                        ++out->element_count;
                        for (size_t i = factor_index; i < a->element_count; ++i)
                        {
                            if (factors[factor_index])
                            {
                                out->elements[out->element_count] =
                                    number_copy(output_stack, factors[factor_index]);
                                ++out->element_count;
                            }
                        }
                        break;
                    }
                }
                ++factor_index;
                if (factor_index == a->element_count)
                {
                    out->elements[out->element_count] =
                        number_copy(output_stack, factors[factor_index]);
                    ++out->element_count;
                    break;
                }
            }
            out->minimal_polynomial = product_minimal_polynomial(output_stack, local_stack, out,
                a->minimal_polynomial, b->minimal_polynomial);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        case '*':
        {
            void*local_stack_savepoint = local_stack->cursor;
            struct Number*out = a;
            for (size_t i = b->element_count; i-- > 1;)
            {
                out = number_multiply(local_stack, output_stack, out, b->elements[i]);
            }
            out = number_multiply(output_stack, local_stack, out, b->elements[0]);
            local_stack->cursor = local_stack_savepoint;
            return out;
        }
        }
        break;
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Number*out = number_rational_initialize(local_stack, &rational_zero);
        for (size_t i = a->element_count; i-- > 1;)
        {
            out = number_add(local_stack, output_stack, out,
                number_multiply(local_stack, output_stack, b, a->elements[i]));
        }
        out = number_add(output_stack, local_stack, out,
            number_multiply(local_stack, output_stack, b, a->elements[0]));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    }
    return number_multiply(output_stack, local_stack, b, a);
}