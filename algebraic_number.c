#include "declarations.h"

int8_t degree_compare(size_t a[2], size_t b[2])
{
    for (size_t i = 0; i < 2; ++i)
    {
        if (a[i] < b[i])
        {
            return -1;
        }
        else if (a[i] > b[i])
        {
            return 1;
        }
    }
    return 0;
}

struct AlgebraicNumber*pool_algebraic_number_allocate(struct Stack*output_stack,
    struct AlgebraicNumber*free_list)
{
    if (free_list)
    {
        struct AlgebraicNumber*out = free_list;
        free_list = free_list->next_term;
        return out;
    }
    return STACK_SLOT_ALLOCATE(output_stack, struct AlgebraicNumber);
}

void algebraic_number_move_coefficients_to_stack(struct Stack*output_stack,
    struct AlgebraicNumber*a)
{
    while (a)
    {
        rational_copy_to_stack(output_stack, a->term_coefficient);
        a = a->next_term;
    }
}

//There must not be anything on output_stack after a. Leaves the coefficients on local_stack.
void algebraic_number_add_to_a_in_place(struct Stack*output_stack, struct Stack*local_stack,
    struct AlgebraicNumber**a, struct AlgebraicNumber*b)
{
    struct AlgebraicNumber*free_list = 0;
    struct AlgebraicNumber*a_term = *a;
    struct AlgebraicNumber*previous_term = 0;
    while (b)
    {
        int8_t degree_comparison = degree_compare(a_term->generator_degrees, b->generator_degrees);
        while (degree_comparison < 0)
        {
            if (a_term->next_term)
            {
                previous_term = a_term;
                a_term = a_term->next_term;
                degree_comparison = degree_compare(a_term->generator_degrees, b->generator_degrees);
            }
            else
            {
                a_term->next_term = pool_algebraic_number_allocate(output_stack, free_list);
                while (b)
                {
                    a_term = a_term->next_term;
                    memcpy(a_term, b, sizeof(struct AlgebraicNumber));
                    b = b->next_term;
                    a_term->next_term = pool_algebraic_number_allocate(output_stack, free_list);
                }
                output_stack->cursor = a_term->next_term;
                a_term->next_term = 0;
                return;
            }
        }
        if (degree_comparison == 0)
        {
            a_term->term_coefficient = rational_add(local_stack, output_stack,
                a_term->term_coefficient, b->term_coefficient);
            if (!a_term->term_coefficient->numerator->value_count)
            {
                if (previous_term)
                {
                    previous_term->next_term = a_term->next_term;
                }
                else
                {
                    *a = a_term->next_term;
                }
                a_term->next_term = free_list;
                free_list = a_term;
            }
        }
        else
        {
            struct AlgebraicNumber*term_to_insert =
                pool_algebraic_number_allocate(output_stack, free_list);
            memcpy(term_to_insert, b, sizeof(struct AlgebraicNumber));
            term_to_insert->next_term = a_term;
            if (previous_term)
            {
                previous_term->next_term = term_to_insert;
            }
            else
            {
                *a = term_to_insert;
            }
        }
        b = b->next_term;
    }
}

struct AlgebraicNumber*algebraic_number_add(struct Stack*output_stack,
    struct Stack*local_stack, struct AlgebraicNumber*a, struct AlgebraicNumber*b)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct AlgebraicNumber*out = STACK_SLOT_ALLOCATE(output_stack, struct AlgebraicNumber);
    struct AlgebraicNumber*out_term = out;
    out_term->next_term = out_term;
    while (a)
    {
        out_term = out_term->next_term;
        memcpy(out_term, a, sizeof(struct AlgebraicNumber));
        a = a->next_term;
        out_term->next_term = STACK_SLOT_ALLOCATE(output_stack, struct AlgebraicNumber);
    }
    output_stack->cursor = out_term->next_term;
    out_term->next_term = 0;    
    algebraic_number_add_to_a_in_place(output_stack, local_stack, &out, b);
    algebraic_number_move_coefficients_to_stack(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

struct AlgebraicNumber*algebraic_number_multiply(struct Stack*output_stack,
    struct Stack*local_stack, struct AlgebraicNumber*a, struct AlgebraicNumber*b,
    struct RationalPolynomial*generator_annulling_polynomials[2])
{
    void*local_stack_savepoint = local_stack->cursor;
    struct AlgebraicNumber*out = 0;
    while (a)
    {
        while (b)
        {
            struct Rational*coefficient = rational_multiply(local_stack, output_stack,
                a->term_coefficient, b->term_coefficient);
            if (coefficient->numerator->value_count)
            {
                struct AlgebraicNumber*product_component =
                    STACK_SLOT_ALLOCATE(local_stack, struct AlgebraicNumber);
                memset(product_component, 0, sizeof(struct AlgebraicNumber));
                product_component->term_coefficient = coefficient;
                for (size_t i = 0; i < 2; ++i)
                {
                    size_t component_generator_degree =
                        a->generator_degrees[i] + b->generator_degrees[i];
                    if (component_generator_degree <
                        generator_annulling_polynomials[i]->coefficient_count - 1)
                    {
                        struct AlgebraicNumber*product_component_term = product_component;
                        while (product_component_term)
                        {
                            product_component_term->generator_degrees[i] =
                                component_generator_degree;
                        }
                    }
                    else
                    {
                        size_t reduced_degree = component_generator_degree -
                            generator_annulling_polynomials[i]->coefficient_count + 1;
                        struct AlgebraicNumber*component_factor =
                            STACK_SLOT_ALLOCATE(local_stack, struct AlgebraicNumber);
                        struct AlgebraicNumber*component_factor_term = component_factor;
                        component_factor_term->next_term = component_factor_term;
                        for (size_t j = 0;
                            j < generator_annulling_polynomials[i]->coefficient_count - 1; ++j)
                        {
                            coefficient = generator_annulling_polynomials[i]->coefficients[j];
                            if (coefficient->numerator->value_count)
                            {
                                component_factor_term = component_factor_term->next_term;
                                memset(component_factor_term, 0, sizeof(struct AlgebraicNumber));
                                component_factor_term->generator_degrees[i] = j + reduced_degree;
                                component_factor_term->term_coefficient =
                                    rational_negative(local_stack, coefficient);
                                component_factor_term->next_term =
                                    STACK_SLOT_ALLOCATE(local_stack, struct AlgebraicNumber);
                            }
                        }
                        component_factor_term->next_term = 0;
                        product_component = algebraic_number_multiply(local_stack, output_stack,
                            product_component, component_factor, generator_annulling_polynomials);
                    }
                }
                algebraic_number_add_to_a_in_place(local_stack, output_stack, &out,
                    product_component);
            }
            b = b->next_term;
        }
        a = a->next_term;
    }
    algebraic_number_move_coefficients_to_stack(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}