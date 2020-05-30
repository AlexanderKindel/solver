#include "declarations.h"

struct Number**binary_sum_or_product_get_conjugates(struct Number*(operation)(struct Stack*restrict,
    struct Stack*restrict, struct Number*, struct Number*), struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number**a_conjugates, size_t a_conjugate_count,
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
            void*output_stack_savepoint = output_stack->cursor;
            out[conjugate_count] = number_eliminate_linear_dependencies(output_stack, local_stack,
                operation(local_stack, output_stack, a_conjugates[i], b_conjugates[j]));
            if (rational_polynomial_equals(minimal_polynomial,
                out[conjugate_count]->minimal_polynomial))
            {
                ++conjugate_count;
                if (conjugate_count == minimal_polynomial->coefficient_count - 1)
                {
                    local_stack->cursor = local_stack_savepoint;
                    return out;
                }
            }
            else
            {
                output_stack->cursor = output_stack_savepoint;
            }
        }
    }
    crash("Not enough conjugates found.");
}

struct Number**sum_or_product_get_conjugates(struct Number*(operation)(struct Stack*restrict,
    struct Stack*restrict, struct Number*, struct Number*),
    struct RationalPolynomial*(get_minimal_polynomial)(struct Stack*restrict, struct Stack*restrict,
        struct Number*, struct RationalPolynomial*, struct RationalPolynomial*),
    struct Stack*restrict output_stack, struct Stack*restrict local_stack, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number**out = number_get_conjugates(local_stack, output_stack, a->elements[0]);
    struct RationalPolynomial*left_minimal_polynomial = a->elements[0]->minimal_polynomial;
    size_t original_element_count = a->element_count;
    a->element_count = 2;
    while (a->element_count < original_element_count)
    {
        struct Number*right_element = a->elements[a->element_count - 1];
        struct RationalPolynomial*new_left_minimal_polynomial =
            get_minimal_polynomial(local_stack, output_stack, a, left_minimal_polynomial,
                right_element->minimal_polynomial);
        out = binary_sum_or_product_get_conjugates(operation, local_stack, output_stack, out,
            left_minimal_polynomial->coefficient_count - 1,
            number_get_conjugates(local_stack, output_stack, right_element),
            right_element->minimal_polynomial->coefficient_count - 1, new_left_minimal_polynomial);
        left_minimal_polynomial = new_left_minimal_polynomial;
        ++a->element_count;
    }
    out = binary_sum_or_product_get_conjugates(operation, output_stack, local_stack, out,
        left_minimal_polynomial->coefficient_count - 1,
        number_get_conjugates(local_stack, output_stack, a->elements[a->element_count - 1]),
        a->elements[a->element_count - 1]->minimal_polynomial->coefficient_count - 1,
        a->minimal_polynomial);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number**number_get_conjugates(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a)
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
        struct Number**radicand_conjugates =
            number_get_conjugates(local_stack, output_stack, a->radicand);
        for (size_t i = 0; i < a->radicand->minimal_polynomial->coefficient_count - 1; ++i)
        {
            radicand_conjugates[i] =
                number_take_root(local_stack, output_stack, radicand_conjugates[i], a->index);
        }
        size_t roots_of_unity_count = integer_to_size_t(a->index);
        struct Number**roots_of_unity = get_roots_of_unity(output_stack, local_stack, a->index);
        struct Number**out =
            binary_sum_or_product_get_conjugates(number_multiply, output_stack, local_stack,
                radicand_conjugates, a->radicand->minimal_polynomial->coefficient_count - 1,
                roots_of_unity, roots_of_unity_count, a->minimal_polynomial);
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    case '*':
        return sum_or_product_get_conjugates(number_multiply, product_get_minimal_polynomial,
            output_stack, local_stack, a);
    case '+':
        return sum_or_product_get_conjugates(number_add, sum_get_minimal_polynomial, output_stack,
            local_stack, a);
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_get_reciprocal(struct Stack*restrict output_stack,
    struct Stack*restrict local_stack, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Number*out;
    switch (a->operation)
    {
    case 'r':
    {
        if (!a->value.numerator->value_count)
        {
            puts("Tried to divide by 0.");
            return 0;
        }
        struct Rational value_reciprocal = rational_get_reciprocal(local_stack, &a->value);
        out = number_rational_initialize(output_stack, &value_reciprocal);
        break;
    }
    case '^':
        out = number_divide(output_stack, local_stack,
            number_take_root(local_stack, output_stack,
                number_integer_exponentiate(local_stack, output_stack, a->radicand,
                    integer_add(local_stack, a->index, INT(1, -1))), a->index), a->radicand);
        break;
    case '*':
    {
        out = &number_one;
        for (size_t i = a->element_count; i-- > 1;)
        {
            out = number_multiply(local_stack, output_stack, out,
                number_get_reciprocal(local_stack, output_stack, a->elements[i]));
        }
        out = number_multiply(output_stack, local_stack, out,
            number_get_reciprocal(local_stack, output_stack, a->elements[0]));
        break;
    }
    case '+':
    {
        a = number_eliminate_linear_dependencies(local_stack, output_stack, a);
        struct Number**conjugates = number_get_conjugates(local_stack, output_stack, a);
        out = &number_one;
        for (size_t i = 1; i < a->minimal_polynomial->coefficient_count - 1; ++i)
        {
            out = number_multiply(local_stack, output_stack, out, conjugates[i]);
        }
        struct Rational coefficient =
            rational_get_reciprocal(local_stack, a->minimal_polynomial->coefficients);
        if (a->minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient.numerator->sign *= -1;
        }
        out = number_rational_multiply(output_stack, local_stack, out, &coefficient);
        break;
    }
    default:
        crash("Number operation not recognized.");
    }
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct Number*number_divide(struct Stack*restrict output_stack, struct Stack*restrict local_stack,
    struct Number*dividend, struct Number*divisor)
{
    struct Number*reciprocal = number_get_reciprocal(local_stack, output_stack, divisor);
    if (!reciprocal)
    {
        return 0;
    }
    return number_multiply(output_stack, local_stack, dividend, reciprocal);
}