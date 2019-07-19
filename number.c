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
            out->field_form->left_term_in_terms_of_generator->coefficients[0]->numerator = &zero;
            out->field_form->left_term_in_terms_of_generator->coefficients[0]->denominator =
                pool_integer_initialize(pool_set, 1, 1);
            out->field_form->left_term_in_terms_of_generator->coefficients[1]->numerator =
                pool_integer_initialize(pool_set, 1, 1);
            out->field_form->left_term_in_terms_of_generator->coefficients[1]->denominator =
                pool_integer_initialize(pool_set, 1, 1);
            out->field_form->right_term_in_terms_of_generator =
                polynomial_allocate(&polynomial_stack, 1);
            out->field_form->right_term_in_terms_of_generator->coefficients[0] =
                POOL_VALUE_ALLOCATE(pool_set, struct Rational);
            out->field_form->right_term_in_terms_of_generator->coefficients[0]->numerator =
                integer_copy_to_pool(pool_set, b->value.numerator);
            out->field_form->right_term_in_terms_of_generator->coefficients[0]->denominator =
                integer_copy_to_pool(pool_set, b->value.denominator);
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
    struct Rational*interval_size = &rational_one;
    while (true)
    {
        for (int i = 0; i < candidate_count;)
        {
            struct Float*real_estimate_min;
            struct Float*real_estimate_max;
            struct Float*imaginary_estimate_min;
            struct Float*imaginary_estimate_max;
            rational_polynomial_estimate_evaluation(pool_set, stack_a, stack_b, &real_estimate_min,
                &real_estimate_max, &imaginary_estimate_min, &imaginary_estimate_max, candidates[i],
                a, interval_size);
            if (float_compare(stack_a, stack_b, &float_zero, real_estimate_min) < 0 ||
                float_compare(stack_a, stack_b, real_estimate_max, &float_zero) < 0 ||
                float_compare(stack_a, stack_b, &float_zero, imaginary_estimate_min) < 0 ||
                float_compare(stack_a, stack_b, imaginary_estimate_max, &float_zero) < 0)
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
            algebraic_form, (struct RationalPolynomial*[2]) {
            a->left->minimal_polynomial,
                a->right->minimal_polynomial
        });
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