#include "declarations.h"

struct Number*number_add(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a, struct Number*b)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->value_count == 0)
        {
            number_free(pool_set, a);
            return b;
        }
        if (b->operation == 'r')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct Number*out = number_allocate(pool_set);
            out->operation = 'r';
            out->value = *rational_add(stack_a, stack_b, &a->value, &b->value);
            rational_move_to_pool(pool_set, &out->value);
            number_free(pool_set, a);
            number_free(pool_set, b);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        break;
    case '^':
    case '*':
        switch (b->operation)
        {
        case 'r':
        {
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->left = a;
            out->right = b;
            out->field_form = POOL_VALUE_ALLOCATE(pool_set, struct SumFieldForm);
            out->field_form->generator = number_copy(pool_set, a);
            out->field_form->left_term_in_terms_of_generator =
                polynomial_allocate(&polynomial_stack, 2);
            out->field_form->left_term_in_terms_of_generator->coefficients[0] =
                POOL_VALUE_ALLOCATE(pool_set, struct Rational);
            out->field_form->left_term_in_terms_of_generator->coefficients[0] =
                rational_copy_to_pool(pool_set, &rational_zero);
            out->field_form->left_term_in_terms_of_generator->coefficients[1] =
                rational_copy_to_pool(pool_set, &rational_one);
            out->field_form->right_term_in_terms_of_generator =
                polynomial_allocate(&polynomial_stack, 1);
            out->field_form->right_term_in_terms_of_generator->coefficients[0] =
                rational_copy_to_pool(pool_set, &b->value);
            return out;
        }
        case '^':
        case '*':
        {
            //Placeholder.
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
        break;
    case '+':
        switch (b->operation)
        {
        case 'r':
        {
            //Placeholder.
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '^':
        {
            //Placeholder.
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '*':
        {
            //Placeholder.
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        case '+':
        {
            //Placeholder.
            struct Number*out = number_allocate(pool_set);
            out->operation = '+';
            out->left = a;
            out->right = b;
            return out;
        }
        }
    }
    return number_add(pool_set, stack_a, stack_b, b, a);
}

struct Number*number_rational_multiply(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*b)
{
    if (b->numerator->value_count == 0)
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value.numerator = &zero;
        out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
        number_free(pool_set, a);
        return out;
    }
    if (integer_equals(b->numerator, &one) && integer_equals(b->denominator, &one))
    {
        return a;
    }
    switch (a->operation)
    {
    case 'r':
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value = *rational_multiply(stack_a, stack_b, &a->value, b);
        rational_move_to_pool(pool_set, &out->value);
        number_free(pool_set, a);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '^':
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = '*';
        out->left = number_allocate(pool_set);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy_to_pool(pool_set, b->numerator);
        out->left->value.denominator = integer_copy_to_pool(pool_set, b->denominator);
        out->right = a;
        return out;
    }
    case '*':
    {
        if (a->left->operation == 'r')
        {
            struct Number*out = number_multiply(pool_set, stack_a, stack_b,
                number_rational_multiply(pool_set, stack_a, stack_b, a->left, b), a->right);
            number_node_free(pool_set, a);
            return out;
        }
        struct Number*out = number_allocate(pool_set);
        out->operation = '*';
        out->left = number_allocate(pool_set);
        out->left->operation = 'r';
        out->left->value.numerator = integer_copy_to_pool(pool_set, b->numerator);
        out->left->value.denominator = integer_copy_to_pool(pool_set, b->denominator);
        out->right = a;
        return out;
    }
    case '+':
    {
        struct Number*out = number_allocate(pool_set);
        out->operation = '+';
        out->left = number_rational_multiply(pool_set, stack_a, stack_b, a->left, b);
        out->right = number_rational_multiply(pool_set, stack_a, stack_b, a->right, b);
        number_node_free(pool_set, a);
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_consolidate_product(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Number*b)
{
    if (b->operation == 'r')
    {
        struct Number*out = number_rational_multiply(pool_set, stack_a, stack_b, a, &b->value);
        number_free(pool_set, b);
        return out;
    }
    switch (a->operation)
    {
    case 'r':
    {
        struct Number*out = number_rational_multiply(pool_set, stack_a, stack_b, b, &a->value);
        number_free(pool_set, a);
        return out;
    }
    case '^':
    {
        if (b->operation == '^')
        {
            void*stack_a_savepoint = stack_a->cursor;
            struct ExtendedGCDInfo info;
            integer_extended_gcd(stack_a, stack_b, &info, a->right->value.denominator,
                b->right->value.denominator);
            struct Integer*product_index = integer_euclidean_quotient(stack_a, stack_b,
                integer_multiply(stack_a, stack_b, a->right->value.denominator,
                    b->right->value.denominator), info.gcd);
            struct Number*radicand_factor_a = number_exponentiate(pool_set, stack_a, stack_b,
                number_copy(pool_set, a->left), &(struct Rational){info.b_over_gcd, &one});
            struct Number*radicand_factor_b = number_exponentiate(pool_set, stack_a, stack_b,
                number_copy(pool_set, b->left), &(struct Rational){info.a_over_gcd, &one});
            struct RationalInterval radicand_factor_a_argument_estimate;
            number_rational_argument_estimate(pool_set, stack_a, stack_b,
                &radicand_factor_a_argument_estimate, radicand_factor_a, &rational_one);
            struct RationalInterval radicand_factor_b_argument_estimate;
            number_rational_argument_estimate(pool_set, stack_a, stack_b,
                &radicand_factor_b_argument_estimate, radicand_factor_b, &rational_one);
            if (!integer_equals(a->right->value.denominator, b->right->value.denominator))
            {
                struct Rational*radicand_factor_argument_interval_size =
                    rational_integer_divide(stack_a, stack_b, &pi_estimate_min, product_index);
                struct RationalInterval a_argument_estimate;
                number_rational_argument_estimate(pool_set, stack_a, stack_b, &a_argument_estimate,
                    a, radicand_factor_argument_interval_size);
                struct RationalInterval b_argument_estimate;
                number_rational_argument_estimate(pool_set, stack_a, stack_b, &b_argument_estimate,
                    b, radicand_factor_argument_interval_size);
                if (rational_compare(stack_a, stack_b,
                    rational_integer_divide(stack_a, stack_b,
                        radicand_factor_a_argument_estimate.max, product_index),
                    a_argument_estimate.min) < 0 ||
                    rational_compare(stack_a, stack_b,
                        rational_integer_divide(stack_a, stack_b,
                            radicand_factor_b_argument_estimate.max, product_index),
                        b_argument_estimate.min) < 0)
                {
                    number_free(pool_set, radicand_factor_a);
                    number_free(pool_set, radicand_factor_b);
                    stack_a->cursor = stack_a_savepoint;
                    return number_product_consolidation_failed;
                }
            }
            struct Number*product_radicand =
                number_multiply(pool_set, stack_a, stack_b, radicand_factor_a, radicand_factor_b);
            struct Number*out = number_exponentiate(pool_set, stack_a, stack_b, product_radicand,
                &(struct Rational){&one, product_index});
            struct RationalInterval product_radicand_rational_argument_estimate;
            number_rational_argument_estimate(pool_set, stack_a, stack_b,
                &product_radicand_rational_argument_estimate, product_radicand, &pi_estimate_min);
            if (rational_compare(stack_a, stack_b, product_radicand_rational_argument_estimate.max,
                rational_add(stack_a, stack_b, radicand_factor_a_argument_estimate.min,
                    radicand_factor_b_argument_estimate.min)) < 0)
            {
                out = number_multiply(pool_set, stack_a, stack_b, out, number_copy(pool_set,
                    get_roots_of_unity(stack_a, stack_b, product_index)[1]));
            }
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        break;
    }
    case '*':
        switch (b->operation)
        {
        case '^':
        {
            struct Number*factor =
                number_consolidate_product(pool_set, stack_a, stack_b, a->left, b);
            if (factor)
            {
                return number_multiply(pool_set, stack_a, stack_b, a->right, factor);
            }
            factor = number_consolidate_product(pool_set, stack_a, stack_b, a->right, b);
            if (factor)
            {
                return number_multiply(pool_set, stack_a, stack_b, a->left, factor);
            }
            return number_product_consolidation_failed;
        }
        case '*':
        {
            struct Number*out = number_multiply(pool_set, stack_a, stack_b, a->left, b);
            if (!out)
            {
                return 0;
            }
            out = number_multiply(pool_set, stack_a, stack_b, out, a->right);
            number_node_free(pool_set, a);
            return out;
        }
        }
        break;
    case '+':
    {
        struct Number*b_copy = number_copy(pool_set, b);
        struct Number*left = number_multiply(pool_set, stack_a, stack_b, a->left, b);
        if (!left)
        {
            return 0;
        }
        struct Number*right = number_multiply(pool_set, stack_a, stack_b, a->right, b_copy);
        if (!right)
        {
            return 0;
        }
        number_node_free(pool_set, a);
        return number_add(pool_set, stack_a, stack_b, left, right);
    }
    }
    return number_multiply(pool_set, stack_a, stack_b, b, a);
}

struct Number*number_multiply(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a, struct Number*b)
{
    struct Number*out = number_consolidate_product(pool_set, stack_a, stack_b, a, b);
    if (out == number_product_consolidation_failed)
    {
        out = number_allocate(pool_set);
        out->operation = '*';
        out->left = a;
        out->right = b;
    }
    return out;
}

struct Number*number_reciprocal(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        if (!a->value.numerator->value_count)
        {
            puts("Tried to divide by 0.");
            return number_divide_by_zero_error;
        }
        POINTER_SWAP(a->value.numerator, a->value.denominator);
        return a;
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        pool_integer_free(pool_set, a->right->value.numerator);
        struct Rational exponent = { integer_add(stack_a, a->right->value.denominator, &INT(1, -)),
            a->right->value.denominator };
        struct Number*radicand_copy = number_copy(pool_set, a->left);
        struct Number*out = number_divide(pool_set, stack_a, stack_b,
            number_exponentiate(pool_set, stack_a, stack_b, a->left, &exponent), radicand_copy);
        number_free(pool_set, a->right);
        number_node_free(pool_set, a);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '*':
    {
        struct Number*out = number_multiply(pool_set, stack_a, stack_b,
            number_reciprocal(pool_set, stack_a, stack_b, a->left),
            number_reciprocal(pool_set, stack_a, stack_b, a->right));
        number_node_free(pool_set, a);
        return out;
    }
    case '+':
    {
        void*stack_a_savepoint = stack_a->cursor;
        number_calculate_conjugates(pool_set, stack_a, stack_b, a);
        struct Number*out = number_allocate(pool_set);
        out->operation = 'r';
        out->value.numerator = pool_integer_initialize(pool_set, 1, 1);
        out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
        struct Number*conjugate = a->next;
        while (conjugate)
        {
            out = number_multiply(pool_set, stack_a, stack_b, out,
                number_copy(pool_set, conjugate));
            conjugate = conjugate->next;
        }
        number_calculate_minimal_polynomial(pool_set, stack_a, stack_b, a);
        struct Rational*coefficient =
            rational_reciprocal(stack_a, stack_b, a->minimal_polynomial->coefficients[0]);
        if (a->minimal_polynomial->coefficient_count % 2 == 0)
        {
            coefficient->numerator->sign *= -1;
        }
        out = number_rational_multiply(pool_set, stack_a, stack_b, out, coefficient);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}

struct Number*number_divide(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct Number*dividend, struct Number*divisor)
{
    struct Number*reciprocal = number_reciprocal(pool_set, stack_a, stack_b, divisor);
    if (reciprocal == number_divide_by_zero_error)
    {
        return number_divide_by_zero_error;
    }
    return number_multiply(pool_set, stack_a, stack_b, dividend, reciprocal);
}

struct Rational*number_rational_factor(struct Stack*output_stack, struct Stack*local_stack,
    struct Number*a)
{
    switch (a->operation)
    {
    case 'r':
        return &a->value;
    case '^':
        return &rational_one;
    case '*':
        return number_rational_factor(output_stack, local_stack, a->left);
    case '+':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct Rational*left_term_rational_factor =
            number_rational_factor(local_stack, output_stack, a->left);
        struct Rational*right_term_rational_factor =
            number_rational_factor(local_stack, output_stack, a->right);
        struct Rational*out = rational_reduced(output_stack, local_stack,
            integer_gcd(local_stack, output_stack, left_term_rational_factor->numerator,
                right_term_rational_factor->numerator),
            integer_lcm(local_stack, output_stack, left_term_rational_factor->denominator,
                right_term_rational_factor->denominator));
        local_stack->cursor = local_stack_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.\n");
    }
}

struct Number*number_exponentiate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*base, struct Rational*exponent)
{
    if (exponent->numerator->sign < 0)
    {
        base = number_reciprocal(pool_set, stack_a, stack_b, base);
        if (base == number_divide_by_zero_error)
        {
            return number_divide_by_zero_error;
        }
        void*stack_a_savepoint = stack_a->cursor;
        struct Number*out = number_exponentiate(pool_set, stack_a, stack_b, base,
            rational_negative(stack_a, exponent));
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    switch (base->operation)
    {
    case 'r':
    {
        if (base->value.numerator->value_count == 0)
        {
            number_free(pool_set, base);
            struct Number*out = number_allocate(pool_set);
            out->operation = 'r';
            out->value.numerator = &zero;
            out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
            return out;
        }
        void*stack_a_savepoint = stack_a->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            rational_move_from_pool(pool_set, stack_a, &base->value);
            base->value =
                *rational_exponentiate(stack_a, stack_b, &base->value, exponent->numerator);
            rational_move_to_pool(pool_set, &base->value);
            stack_a->cursor = stack_a_savepoint;
            return base;
        }
        if (!integer_equals(base->value.denominator, &one))
        {
            struct Number*new_numerator_base = number_allocate(pool_set);
            new_numerator_base->operation = 'r';
            new_numerator_base->value.denominator = pool_integer_initialize(pool_set, 1, 1);
            new_numerator_base->value.numerator =
                integer_multiply(stack_a, stack_b, base->value.numerator,
                    integer_exponentiate(stack_a, stack_b, base->value.denominator,
                        integer_add(stack_a, exponent->denominator, &INT(1, -))));
            integer_move_to_pool(pool_set, &new_numerator_base->value.numerator);
            struct Number*new_numerator =
                number_exponentiate(pool_set, stack_a, stack_b, new_numerator_base, exponent);
            struct Number*out = number_rational_multiply(pool_set, stack_a, stack_b, new_numerator,
                &(struct Rational){&one, base->value.denominator});
            number_free(pool_set, base);
            stack_a->cursor = stack_a_savepoint;
            return out;
        }
        if (!integer_equals(exponent->denominator, &one))
        {
            integer_move_from_pool(pool_set, stack_a, &base->value.numerator);
            int8_t base_sign = base->value.numerator->sign;
            base->value.numerator->sign = 1;
            struct Factor*factors;
            size_t factor_count = integer_factor(stack_a, stack_b, &factors, base->value.numerator);
            struct Integer*coefficient = &one;
            struct Integer*multiplicity_gcd = exponent->denominator;
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct IntegerDivision multiplicity_reduction;
                integer_euclidean_divide(stack_a, stack_b, &multiplicity_reduction,
                    factors[factor_index].multiplicity, exponent->denominator);
                factors[factor_index].multiplicity = multiplicity_reduction.remainder;
                coefficient = integer_multiply(stack_a, stack_b, coefficient,
                    integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                        multiplicity_reduction.quotient));
                multiplicity_gcd = integer_gcd(stack_a, stack_b, multiplicity_gcd,
                    factors[factor_index].multiplicity);
            }
            if (base_sign > 0)
            {
                struct Integer*reduced_degree = integer_euclidean_quotient(stack_a, stack_b,
                    exponent->denominator, multiplicity_gcd);
                if (integer_equals(reduced_degree, &one))
                {
                    struct Number*out = number_allocate(pool_set);
                    out->operation = 'r';
                    out->value.denominator = pool_integer_initialize(pool_set, 1, 1);
                    out->value.numerator = coefficient;
                    integer_move_to_pool(pool_set, &out->value.numerator);
                    stack_a->cursor = stack_a_savepoint;
                    return out;
                }
                exponent = STACK_SLOT_ALLOCATE(stack_a, struct Rational);
                exponent->numerator = &one;
                exponent->denominator = reduced_degree;
            }
            base->value.numerator = stack_integer_initialize(stack_a, 1, base_sign);
            for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
            {
                struct Integer*reduced_multiplicity = integer_euclidean_quotient(stack_a, stack_b,
                    factors[factor_index].multiplicity, multiplicity_gcd);
                struct Integer*exponentiation = integer_exponentiate(stack_a, stack_b,
                    factors[factor_index].value, reduced_multiplicity);
                base->value.numerator =
                    integer_multiply(stack_a, stack_b, base->value.numerator, exponentiation);
            }
            integer_move_to_pool(pool_set, &base->value.numerator);
            struct Number*number_exponent = number_allocate(pool_set);
            number_exponent->operation = 'r';
            number_exponent->value.numerator = pool_integer_initialize(pool_set, 1, 1);
            number_exponent->value.denominator =
                integer_copy_to_pool(pool_set, exponent->denominator);
            struct Number*surd = number_allocate(pool_set);
            surd->operation = '^';
            surd->left = base;
            surd->right = number_exponent;
            stack_a->cursor = stack_a_savepoint;
            return number_rational_multiply(pool_set, stack_a, stack_b, surd,
                &(struct Rational){coefficient, &one});
        }
        return base;
    }
    case '^':
    {
        void*stack_a_savepoint = stack_a->cursor;
        exponent = rational_multiply(stack_a, stack_b, exponent, &base->right->value);
        number_free(pool_set, base->right);
        struct Number*out = number_exponentiate(pool_set, stack_a, stack_b, base->left, exponent);
        number_node_free(pool_set, base);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    case '*':
    {
        struct Number*out = number_multiply(pool_set, stack_a, stack_b,
            number_exponentiate(pool_set, stack_a, stack_b, base->left, exponent),
            number_exponentiate(pool_set, stack_a, stack_b, base->right, exponent));
        number_node_free(pool_set, base);
        return out;
    }
    case '+':
    {
        void*stack_a_savepoint = stack_a->cursor;
        if (!integer_equals(exponent->numerator, &one))
        {
            struct Number*exponentiation_by_numerator = number_allocate(pool_set);
            exponentiation_by_numerator->operation = 'r';
            exponentiation_by_numerator->value.numerator = pool_integer_initialize(pool_set, 1, 1);
            exponentiation_by_numerator->value.denominator =
                pool_integer_initialize(pool_set, 1, 1);
            struct Integer*numerator = exponent->numerator;
            while (numerator->sign > 0)
            {
                if (numerator->value[0] & 1)
                {
                    exponentiation_by_numerator = number_multiply(pool_set, stack_a, stack_b,
                        exponentiation_by_numerator, number_copy(pool_set, base));
                }
                base =
                    number_multiply(pool_set, stack_a, stack_b, base, number_copy(pool_set, base));
                numerator = integer_half(stack_a, numerator);
            }
            number_free(pool_set, base);
            stack_a->cursor = stack_a_savepoint;
            return number_exponentiate(pool_set, stack_a, stack_b, exponentiation_by_numerator,
                &(struct Rational){&one, exponent->denominator});
        }
        struct Rational*base_rational_factor = number_rational_factor(stack_a, stack_b, base);
        struct Factor*factors;
        size_t factor_count =
            integer_factor(stack_a, stack_b, &factors, base_rational_factor->numerator);
        struct Rational base_cancelling_rational_factor =
            { integer_exponentiate(stack_a, stack_b, base_rational_factor->denominator,
                exponent->denominator), &one };
        struct Rational product_rational_factor = { &one, base_rational_factor->denominator };
        for (size_t factor_index = 0; factor_index < factor_count; ++factor_index)
        {
            struct Integer*reduced_multiplicity = integer_euclidean_quotient(stack_a, stack_b,
                factors[factor_index].multiplicity, exponent->denominator);
            base_cancelling_rational_factor.denominator = integer_multiply(stack_a, stack_b,
                integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                    integer_multiply(stack_a, stack_b, reduced_multiplicity,
                        exponent->denominator)), base_cancelling_rational_factor.denominator);
            product_rational_factor.numerator =
                integer_multiply(stack_a, stack_b, product_rational_factor.numerator,
                    integer_exponentiate(stack_a, stack_b, factors[factor_index].value,
                        reduced_multiplicity));
        }
        struct Number*surd = number_allocate(pool_set);
        surd->operation = '^';
        surd->left = number_rational_multiply(pool_set, stack_a, stack_b, base,
            &base_cancelling_rational_factor);
        surd->right = number_allocate(pool_set);
        surd->right->operation = 'r';
        surd->right->value.numerator = pool_integer_initialize(pool_set, 1, 1);
        surd->right->value.denominator = integer_copy_to_pool(pool_set, exponent->denominator);
        struct Number*out =
            number_rational_multiply(pool_set, stack_a, stack_b, surd, &product_rational_factor);
        stack_a->cursor = stack_a_savepoint;
        return out;
    }
    default:
        crash("Number operation not recognized.");
    }
}