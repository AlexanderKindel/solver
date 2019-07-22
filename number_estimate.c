#include "declarations.h"

struct Rational*rational_place_value(struct Stack*output_stack, struct Stack*local_stack,
    struct Rational*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = &rational_one;
    while (rational_compare(output_stack, local_stack, out, a) > 0)
    {
        out->denominator = integer_doubled(local_stack, out->denominator);
    }
    out->denominator = integer_doubled(output_stack, out->denominator);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_float_estimate_from_rational(rational_estimate_getter get_rational_estimate,
    struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b, struct FloatInterval*out,
    struct Number*a, struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*rational_estimate_interval_size =
        rational_place_value(stack_a, stack_b, interval_size);
    struct RationalInterval rational_estimate;
    get_rational_estimate(pool_set, stack_a, stack_b, &rational_estimate, a,
        rational_estimate_interval_size);
    rational_interval_to_float_interval(stack_a, stack_b, out, &rational_estimate,
        rational_estimate_interval_size);
    float_interval_move_value_to_pool(pool_set, out);
    stack_a->cursor = stack_a_savepoint;
}

bool estimate_needs_refinement(struct PoolSet*pool_set, struct Stack*stack_a, struct Stack*stack_b,
    struct FloatInterval**estimate, struct Rational*interval_size)
{
    if (*estimate)
    {
        void*stack_a_savepoint = stack_a->cursor;
        struct Float*current_interval_size =
            float_subtract(stack_a, stack_b, (*estimate)->max, (*estimate)->min);
        struct Rational*rational_current_interval_size =
            float_to_rational(stack_a, stack_b, current_interval_size);
        if (rational_compare(stack_a, stack_b, interval_size, rational_current_interval_size) >= 0)
        {
            stack_a->cursor = stack_a_savepoint;
            return false;
        }
        else
        {
            float_free(pool_set, (*estimate)->min);
            float_free(pool_set, (*estimate)->max);
        }
        stack_a->cursor = stack_a_savepoint;
    }
    else
    {
        *estimate = pool_value_allocate(pool_set, sizeof(struct FloatInterval));
    }
    return true;
}

struct Rational*number_estimate_max_magnitude(rational_estimate_getter get_estimate,
    struct PoolSet*pool_set, struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval estimate;
    get_estimate(pool_set, local_stack, output_stack, &estimate, a, &rational_one);
    struct Rational*out = rational_max(output_stack, local_stack,
        rational_magnitude(local_stack, estimate.min),
        rational_magnitude(local_stack, estimate.max));
    out = rational_copy_to_stack(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_argument_cosine_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(pool_set, local_stack, output_stack, &argument_estimate, a,
        argument_estimate_interval_size);
    struct RationalInterval argument_min_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_min_cosine, argument_estimate.min,
        argument_estimate_interval_size);
    struct RationalInterval argument_max_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_max_cosine, argument_estimate.max,
        argument_estimate_interval_size);
    out->max = rational_copy_to_stack(output_stack,
        rational_max(output_stack, local_stack, argument_min_cosine.max, argument_max_cosine.max));
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack, argument_estimate.min);
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack, argument_estimate.max);
    struct Rational*pi_estimate_max =
        rational_add(local_stack, output_stack, &pi_estimate_min, &pi_interval_size);
    if (rational_compare(output_stack, local_stack, argument_estimate.min, &pi_estimate_min) <= 0 &&
        rational_compare(output_stack, local_stack, pi_estimate_max, argument_estimate.max) <= 0)
    {
        out->min->numerator = stack_integer_initialize(output_stack, 1, -1);
        out->min->denominator = &one;
    }
    else
    {
        out->min = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
            argument_min_cosine.min, argument_max_cosine.min));
        out->min = rational_copy_to_stack(output_stack, out->min);
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_argument_sine_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(pool_set, local_stack, output_stack, &argument_estimate, a,
        argument_estimate_interval_size);
    struct RationalInterval argument_min_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_min_sine, argument_estimate.min,
        argument_estimate_interval_size);
    struct RationalInterval argument_max_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_max_sine, argument_estimate.max,
        argument_estimate_interval_size);
    struct RationalInterval argument_estimate_multiple =
        { rational_doubled(local_stack, output_stack, argument_estimate.min),
    rational_doubled(local_stack, output_stack, argument_estimate.max) };
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
        argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
        argument_estimate_multiple.max);
    struct Rational*pi_estimate_max =
        rational_add(local_stack, output_stack, &pi_estimate_min, &pi_interval_size);
    if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
        &pi_estimate_min) <= 0 && rational_compare(output_stack, local_stack, pi_estimate_max,
            argument_estimate_multiple.max) <= 0)
    {
        out->max->numerator = stack_integer_initialize(output_stack, 1, 1);
        out->max->denominator = stack_integer_initialize(output_stack, 1, 1);
        out->min = rational_copy_to_stack(output_stack,
            rational_max(output_stack, local_stack, argument_min_sine.min, argument_max_sine.min));
    }
    else
    {
        argument_estimate_multiple.min = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.min, &INT(3, +));
        argument_estimate_multiple.max = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.max, &INT(3, +));
        pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
            argument_estimate_multiple.min);
        pi_shrink_interval_to_one_side_of_value(output_stack, local_stack,
            argument_estimate_multiple.max);
        pi_estimate_max =
            rational_add(local_stack, output_stack, &pi_estimate_min, &pi_interval_size);
        if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
            &pi_estimate_min) <= 0 && rational_compare(output_stack, local_stack, pi_estimate_max,
                argument_estimate_multiple.max) <= 0)
        {
            out->min->numerator = stack_integer_initialize(output_stack, 1, -1);
            out->min->denominator = stack_integer_initialize(output_stack, 1, 1);
            out->max = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
        else
        {
            out->min = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.min, argument_max_sine.min));
            out->max = rational_copy_to_stack(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_rectangular_part_from_polar_form(void(trig_function)(struct PoolSet*, struct Stack*,
    struct Stack*, struct RationalInterval*, struct Number*, struct Rational*),
    struct PoolSet*pool_set, struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval magnitude_estimate;
    number_rational_magnitude_estimate(pool_set, local_stack, output_stack, &magnitude_estimate, a,
        &rational_one);
    struct Rational*factor_interval_size = rational_divide(local_stack, output_stack, interval_size,
        rational_integer_add(local_stack, output_stack, magnitude_estimate.max, &INT(2, +)));
    number_rational_magnitude_estimate(pool_set, local_stack, output_stack, &magnitude_estimate, a,
        factor_interval_size);
    struct RationalInterval trig_value;
    trig_function(pool_set, local_stack, output_stack, &trig_value, a, factor_interval_size);
    if (trig_value.min->numerator->sign >= 0)
    {
        out->min =
            rational_multiply(output_stack, local_stack, trig_value.min, magnitude_estimate.min);
        out->max =
            rational_multiply(output_stack, local_stack, trig_value.max, magnitude_estimate.max);
    }
    else if (trig_value.max->numerator->sign <= 0)
    {
        out->min =
            rational_multiply(output_stack, local_stack, trig_value.min, magnitude_estimate.max);
        out->max =
            rational_multiply(output_stack, local_stack, trig_value.max, magnitude_estimate.min);
    }
    else
    {
        out->min =
            rational_multiply(output_stack, local_stack, trig_value.min, magnitude_estimate.max);
        out->max =
            rational_multiply(output_stack, local_stack, trig_value.max, magnitude_estimate.max);
    }
    local_stack->cursor = local_stack_savepoint;
}

struct Rational*factor_estimate_interval_size(rational_estimate_getter get_factor_a_estimate,
    rational_estimate_getter get_factor_b_estimate, struct PoolSet*pool_set,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*factor_a,
    struct Number*factor_b, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_divide(output_stack, local_stack, interval_size,
        rational_add(local_stack, output_stack,
            number_estimate_max_magnitude(get_factor_a_estimate, pool_set, local_stack,
                output_stack, factor_a),
            rational_integer_add(local_stack, output_stack,
                number_estimate_max_magnitude(get_factor_b_estimate, pool_set, local_stack,
                    output_stack, factor_b), &one)));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_estimate_part_sum(float_estimate_getter get_term_estimate, struct PoolSet*pool_set,
    struct Stack*stack_a, struct Stack*stack_b, struct FloatInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    struct Rational*term_estimate_interval_size =
        rational_integer_divide(stack_a, stack_b, interval_size, &INT(2, +));
    struct FloatInterval*left_term_estimate =
        get_term_estimate(pool_set, stack_a, stack_b, a->left, term_estimate_interval_size);
    struct FloatInterval*right_term_estimate =
        get_term_estimate(pool_set, stack_a, stack_b, a->right, term_estimate_interval_size);
    out->min = float_add(stack_a, stack_b, left_term_estimate->min, right_term_estimate->min);
    out->max = float_add(stack_a, stack_b, left_term_estimate->max, right_term_estimate->max);
    float_interval_move_value_to_pool(pool_set, out);
    stack_a->cursor = stack_a_savepoint;
}

//When a takes the union variant with the real_part_estimate field, the return value is that field,
//allocated in pool_set. Otherwise, the return value is allocated on output_stack.
struct FloatInterval*number_float_real_part_estimate(struct PoolSet*pool_set,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == 'r')
    {
        struct FloatInterval*magnitude_estimate =
            number_float_magnitude_estimate(pool_set, output_stack, local_stack, a, interval_size);
        if (a->value.numerator->sign < 0)
        {
            struct FloatInterval*real_part_estimate =
                STACK_SLOT_ALLOCATE(output_stack, struct FloatInterval);
            real_part_estimate->min = float_negative(output_stack, magnitude_estimate->max);
            real_part_estimate->max = float_negative(output_stack, magnitude_estimate->min);
            return real_part_estimate;
        }
        else
        {
            return magnitude_estimate;
        }
    }
    if (!estimate_needs_refinement(pool_set, output_stack, local_stack, &a->real_part_estimate,
        interval_size))
    {
        return a->real_part_estimate;
    }
    switch (a->operation)
    {
    case '^':
        number_float_estimate_from_rational(number_rational_real_part_estimate, pool_set,
            output_stack, local_stack, a->real_part_estimate, a, interval_size);
        return a->real_part_estimate;
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        interval_size = rational_half(local_stack, output_stack, interval_size);
        struct Rational*real_part_interval_size =
            factor_estimate_interval_size(number_rational_real_part_estimate,
                number_rational_real_part_estimate, pool_set, local_stack, output_stack, a->left,
                a->right, interval_size);
        struct FloatInterval real_part_product;
        float_interval_multiply(local_stack, output_stack, &real_part_product,
            number_float_real_part_estimate(pool_set, local_stack, output_stack, a->left,
                real_part_interval_size),
            number_float_real_part_estimate(pool_set, local_stack, output_stack, a->right,
                real_part_interval_size));
        struct Rational*imaginary_part_interval_size =
            factor_estimate_interval_size(number_rational_imaginary_part_estimate,
                number_rational_imaginary_part_estimate, pool_set, local_stack, output_stack,
                a->left, a->right, interval_size);
        struct FloatInterval imaginary_part_product;
        float_interval_multiply(local_stack, output_stack, &imaginary_part_product,
            number_float_imaginary_part_estimate(pool_set, local_stack, output_stack, a->left,
                imaginary_part_interval_size),
            number_float_real_part_estimate(pool_set, local_stack, output_stack, a->right,
                imaginary_part_interval_size));
        a->real_part_estimate->min = float_subtract(local_stack, output_stack,
            real_part_product.min, imaginary_part_product.max);
        a->real_part_estimate->max = float_subtract(local_stack, output_stack,
            real_part_product.max, imaginary_part_product.min);
        float_interval_move_value_to_pool(pool_set, a->real_part_estimate);
        local_stack->cursor = local_stack_savepoint;
        return a->real_part_estimate;
    }
    case '+':
        number_estimate_part_sum(number_float_real_part_estimate, pool_set, output_stack,
            local_stack, a->real_part_estimate, a, interval_size);
        return a->real_part_estimate;
    default:
        crash("Number operation not recognized.\n");
    }
}

void number_rational_real_part_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_cosine_estimate, pool_set,
            output_stack, local_stack, out, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval*float_estimate =
            number_float_real_part_estimate(pool_set, local_stack, output_stack, a, interval_size);
        out->min = float_to_rational(output_stack, local_stack, float_estimate->min);
        out->max = float_to_rational(output_stack, local_stack, float_estimate->max);
        local_stack->cursor = local_stack_savepoint;
    }
}

//When a takes the union variant with the imaginary_part_estimate field, the return value is that
//field, allocated in pool_set. Otherwise, the return value is allocated on output_stack.
struct FloatInterval*number_float_imaginary_part_estimate(struct PoolSet*pool_set,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == 'r')
    {
        struct FloatInterval*imaginary_part_estimate =
            STACK_SLOT_ALLOCATE(output_stack, struct FloatInterval);
        imaginary_part_estimate->min = STACK_SLOT_ALLOCATE(output_stack, struct Float);
        imaginary_part_estimate->max = STACK_SLOT_ALLOCATE(output_stack, struct Float);
        imaginary_part_estimate->min->significand = &zero;
        imaginary_part_estimate->min->exponent = &zero;
        imaginary_part_estimate->max->significand = &zero;
        imaginary_part_estimate->max->exponent = &zero;
        return imaginary_part_estimate;
    }
    if (!estimate_needs_refinement(pool_set, output_stack, local_stack, &a->imaginary_part_estimate,
        interval_size))
    {
        return a->imaginary_part_estimate;
    }
    switch (a->operation)
    {
    case '^':
    {
        number_float_estimate_from_rational(number_rational_imaginary_part_estimate, pool_set,
            output_stack, local_stack, a->imaginary_part_estimate, a, interval_size);
        return a->imaginary_part_estimate;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        interval_size = rational_half(local_stack, output_stack, interval_size);
        struct Rational*term_a_interval_size =
            factor_estimate_interval_size(number_rational_real_part_estimate,
                number_rational_imaginary_part_estimate, pool_set, local_stack, output_stack,
                a->left, a->right, interval_size);
        struct FloatInterval term_a;
        float_interval_multiply(local_stack, output_stack, &term_a,
            number_float_real_part_estimate(pool_set, local_stack, output_stack, a->left,
                term_a_interval_size),
            number_float_imaginary_part_estimate(pool_set, local_stack, output_stack, a->right,
                term_a_interval_size));
        struct Rational*term_b_interval_size =
            factor_estimate_interval_size(number_rational_imaginary_part_estimate,
                number_rational_real_part_estimate, pool_set, local_stack, output_stack,
                a->left, a->right, interval_size);
        struct FloatInterval term_b;
        float_interval_multiply(local_stack, output_stack, &term_b,
            number_float_imaginary_part_estimate(pool_set, local_stack, output_stack, a->left,
                term_b_interval_size),
            number_float_real_part_estimate(pool_set, local_stack, output_stack, a->right,
                term_b_interval_size));
        a->imaginary_part_estimate->min =
            float_add(local_stack, output_stack, term_a.min, term_b.min);
        a->imaginary_part_estimate->max =
            float_add(local_stack, output_stack, term_a.max, term_b.max);
        float_interval_move_value_to_pool(pool_set, a->imaginary_part_estimate);
        local_stack->cursor = local_stack_savepoint;
        return a->imaginary_part_estimate;
    }
    case '+':
        number_estimate_part_sum(number_float_imaginary_part_estimate, pool_set, output_stack,
            local_stack, a->imaginary_part_estimate, a, interval_size);
        return a->imaginary_part_estimate;
    default:
        crash("Number operation not recognized.\n");
    }
}

void number_rational_imaginary_part_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_sine_estimate, pool_set,
            output_stack, local_stack, out, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval*float_estimate = number_float_imaginary_part_estimate(pool_set,
            local_stack, output_stack, a, interval_size);
        out->min = float_to_rational(output_stack, local_stack, float_estimate->min);
        out->max = float_to_rational(output_stack, local_stack, float_estimate->max);
        local_stack->cursor = local_stack_savepoint;
    }
}

struct FloatInterval*number_rectangular_part_estimate_for_magnitude_estimate(
    float_estimate_getter get_part_float_estimate,
    rational_estimate_getter get_part_rational_estimate, struct PoolSet*pool_set,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval*out = get_part_float_estimate(pool_set, output_stack, local_stack, a,
        rational_divide(local_stack, output_stack, interval_size,
            rational_integer_add(local_stack, output_stack,
                rational_integer_multiply(local_stack, output_stack,
                    number_estimate_max_magnitude(get_part_rational_estimate, pool_set, local_stack,
                        output_stack, a), &INT(4, +)), &INT(2, +))));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_magnitude_estimate_bound_interval(struct Stack*output_stack, struct Stack*local_stack,
    struct Float**out_min, struct Float**out_max, struct Float*real_part_bound,
    struct Float*imaginary_part_bound, struct Rational*interval_size)
{
    float_estimate_root(output_stack, local_stack, out_min, out_max,
        float_add(local_stack, output_stack,
            float_multiply(local_stack, output_stack, real_part_bound, real_part_bound),
            float_multiply(local_stack, output_stack, imaginary_part_bound, imaginary_part_bound)),
        interval_size, &INT(2, +));
}

struct FloatInterval*number_float_magnitude_estimate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*interval_size)
{
    void*stack_a_savepoint = stack_a->cursor;
    if (a->operation == 'r')
    {
        if (!a->magnitude_estimate)
        {
            a->magnitude_estimate = pool_value_allocate(pool_set, sizeof(struct FloatInterval));
            struct IntegerDivision division;
            integer_euclidean_divide(stack_a, stack_b, &division,
                integer_magnitude(stack_a, a->value.numerator), a->value.denominator);
            a->magnitude_estimate_denominator = &one;
            a->magnitude_estimate->min->significand = division.quotient;
            a->magnitude_estimate->min->exponent = &zero;
            a->magnitude_estimate_remainder = division.remainder;
        }
        else
        {
            float_move_from_pool(pool_set, stack_a, &a->magnitude_estimate->min);
            float_move_from_pool(pool_set, stack_a, &a->magnitude_estimate->max);
            integer_move_from_pool(pool_set, stack_a, &a->magnitude_estimate_denominator);
            integer_move_from_pool(pool_set, stack_a, &a->magnitude_estimate_remainder);
        }
        rational_continue_float_estimate(stack_a, stack_b, &a->argument_estimate->min,
            &a->argument_estimate->max, &a->magnitude_estimate_denominator,
            &a->magnitude_estimate_remainder, &a->value, interval_size);
        float_interval_move_value_to_pool(pool_set, a->magnitude_estimate);
        integer_move_to_pool(pool_set, &a->magnitude_estimate_denominator);
        integer_move_to_pool(pool_set, &a->magnitude_estimate_remainder);
        stack_a->cursor = stack_a_savepoint;
        return a->magnitude_estimate;
    }
    if (!estimate_needs_refinement(pool_set, stack_a, stack_b, &a->magnitude_estimate,
        interval_size))
    {
        return a->magnitude_estimate;
    }
    switch (a->operation)
    {
    case '^':
    {
        struct FloatInterval*radicand_magnitude_estimate =
            number_float_magnitude_estimate(pool_set, stack_a, stack_b, a->left,
                rational_integer_divide(stack_a, stack_b,
                    rational_exponentiate(stack_a, stack_b, interval_size,
                        a->right->value.denominator),
                    &INT(3, +)));
        struct Rational*bound_interval_size =
            rational_integer_divide(stack_a, stack_b, interval_size, &INT(3, +));
        struct Float*unused_estimate_bound;
        float_estimate_root(stack_a, stack_b, &a->magnitude_estimate->min, &unused_estimate_bound,
            radicand_magnitude_estimate->min, bound_interval_size, a->right->value.denominator);
        float_estimate_root(stack_a, stack_b, &unused_estimate_bound, &a->magnitude_estimate->max,
            radicand_magnitude_estimate->max, bound_interval_size, a->right->value.denominator);
        float_interval_move_value_to_pool(pool_set, a->magnitude_estimate);
        stack_a->cursor = stack_a_savepoint;
        return a->magnitude_estimate;
    }
    case '*':
    {
        struct RationalInterval left_factor_magnitude_estimate;
        number_rational_magnitude_estimate(pool_set, stack_a, stack_b,
            &left_factor_magnitude_estimate, a->left, &rational_one);
        struct RationalInterval right_factor_magnitude_estimate;
        number_rational_magnitude_estimate(pool_set, stack_a, stack_b,
            &right_factor_magnitude_estimate, a->right, &rational_one);
        interval_size = rational_divide(stack_a, stack_b, interval_size,
            rational_multiply(stack_a, stack_b,
                rational_add(stack_a, stack_b, left_factor_magnitude_estimate.max, &rational_one),
                rational_add(stack_a, stack_b, right_factor_magnitude_estimate.max,
                    &rational_one)));
        number_float_magnitude_estimate(pool_set, stack_a, stack_b, a->left, interval_size);
        number_float_magnitude_estimate(pool_set, stack_a, stack_b, a->right, interval_size);
        a->argument_estimate->min = float_multiply(stack_a, stack_b,
            a->left->magnitude_estimate->min, a->right->magnitude_estimate->min);
        a->argument_estimate->max = float_multiply(stack_a, stack_b,
            a->left->magnitude_estimate->max, a->right->magnitude_estimate->max);
        float_interval_move_value_to_pool(pool_set, a->magnitude_estimate);
        stack_a->cursor = stack_a_savepoint;
        return a->magnitude_estimate;
    }
    case '+':
    {
        interval_size = rational_integer_divide(stack_a, stack_b, interval_size, &INT(3, +));
        struct FloatInterval*real_part_estimate =
            number_rectangular_part_estimate_for_magnitude_estimate(number_float_real_part_estimate,
                number_rational_real_part_estimate, pool_set, stack_a, stack_b, a, interval_size);
        struct FloatInterval*imaginary_part_estimate =
            number_rectangular_part_estimate_for_magnitude_estimate(
                number_float_imaginary_part_estimate, number_rational_imaginary_part_estimate,
                pool_set, stack_a, stack_b, a, interval_size);
        struct Float*unused_estimate_bound;
        number_magnitude_estimate_bound_interval(stack_a, stack_b, &a->argument_estimate->min,
            &unused_estimate_bound, real_part_estimate->min, imaginary_part_estimate->min,
            interval_size);
        number_magnitude_estimate_bound_interval(stack_a, stack_b, &unused_estimate_bound,
            &a->argument_estimate->max, real_part_estimate->max, imaginary_part_estimate->max,
            interval_size);
        float_interval_move_value_to_pool(pool_set, a->magnitude_estimate);
        stack_a->cursor = stack_a_savepoint;
        return a->magnitude_estimate;
    }
    default:
        crash("Number operation not recognized.");
    }
}

void number_rational_magnitude_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    float_interval_to_rational_interval(output_stack, local_stack, out,
        number_float_magnitude_estimate(pool_set, output_stack, local_stack, a, interval_size));
}

struct FloatInterval*float_argument_estimate_for_part_sum(struct PoolSet*pool_set,
    struct Stack*stack_a, struct Stack*stack_b, struct Number*a, struct Rational*interval_size)
{
    return number_float_argument_estimate(pool_set, stack_a, stack_b, a, interval_size);
}

struct FloatInterval*number_float_argument_estimate(struct PoolSet*pool_set, struct Stack*stack_a,
    struct Stack*stack_b, struct Number*a, struct Rational*interval_size)
{
    if (!estimate_needs_refinement(pool_set, stack_a, stack_b, &a->argument_estimate,
        interval_size))
    {
        return a->magnitude_estimate;
    }
    if (a->operation == '*')
    {
        number_estimate_part_sum(float_argument_estimate_for_part_sum, pool_set, stack_a, stack_b,
            a->argument_estimate, a, interval_size);
        return a->argument_estimate;
    }
    else
    {
        number_float_estimate_from_rational(number_rational_argument_estimate, pool_set, stack_a,
            stack_b, a->argument_estimate, a, interval_size);
        return a->argument_estimate;
    }
}

void number_rational_argument_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RationalInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->sign < 0)
        {
            pi_estimate(output_stack, local_stack, interval_size);
            out->min = &pi_estimate_min;
            out->max = rational_add(output_stack, local_stack, &pi_estimate_min, &pi_interval_size);
        }
        else
        {
            out->min = STACK_SLOT_ALLOCATE(output_stack, struct Rational);
            out->min->numerator = &zero;
            out->min->denominator = &one;
            out->max = STACK_SLOT_ALLOCATE(output_stack, struct Rational);
            out->max->numerator = &zero;
            out->max->denominator = &one;
        }
        return;
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval radicand_rational_argument_estimate;
        number_rational_argument_estimate(pool_set, local_stack, output_stack,
            &radicand_rational_argument_estimate, a->left,
            rational_integer_multiply(local_stack, output_stack, interval_size,
                a->right->value.denominator));
        out->min = rational_integer_divide(output_stack, local_stack,
            radicand_rational_argument_estimate.min, a->right->value.denominator);
        out->max = rational_integer_divide(output_stack, local_stack,
            radicand_rational_argument_estimate.max, a->right->value.denominator);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        float_interval_to_rational_interval(output_stack, local_stack, out,
            number_float_argument_estimate(pool_set, output_stack, local_stack, a, interval_size));
        local_stack->cursor = local_stack_savepoint;
    }
    case '+':
        crash("number_rational_argument_estimate case not yet implemented.\n");
    }
}

void number_rectangular_estimate(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct Number*a,
    struct Rational*interval_size)
{
    out->real_part_estimate =
        number_float_real_part_estimate(pool_set, output_stack, local_stack, a, interval_size);
    out->imaginary_part_estimate =
        number_float_imaginary_part_estimate(pool_set, output_stack, local_stack, a, interval_size);
}

bool rectangular_estimates_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct RectangularEstimate*a, struct RectangularEstimate*b)
{
    return float_intervals_are_disjoint(stack_a, stack_b, a->real_part_estimate,
            b->real_part_estimate) ||
        float_intervals_are_disjoint(stack_a, stack_b, a->imaginary_part_estimate,
            b->imaginary_part_estimate);
}

void rational_polynomial_estimate_evaluation(struct PoolSet*pool_set, struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct RationalPolynomial*a,
    struct Number*argument, struct Rational*interval_size)
{
    if (a->coefficient_count == 0)
    {
        out->real_part_estimate->min = &float_zero;
        out->real_part_estimate->max = &float_zero;
        out->imaginary_part_estimate->min = &float_zero;
        out->imaginary_part_estimate->max = &float_zero;
        return;
    }
    if (a->coefficient_count == 1)
    {
        rational_float_estimate(output_stack, local_stack, &out->real_part_estimate->min,
            &out->real_part_estimate->max, a->coefficients[0], interval_size);
        out->imaginary_part_estimate->min = &float_zero;
        out->imaginary_part_estimate->max = &float_zero;
        return;
    }
    void*local_stack_savepoint = local_stack->cursor;
    struct Integer*evaluation_real_term_count;
    struct Integer*evaluation_imaginary_term_count;
    struct IntegerDivision division;
    integer_euclidean_divide(local_stack, output_stack, &division,
        integer_from_size_t(local_stack, a->coefficient_count), &INT(2, +));
    if (!division.remainder->value_count)
    {
        evaluation_imaginary_term_count =
            integer_multiply(local_stack, output_stack, division.quotient, division.quotient);
        evaluation_real_term_count = integer_add(local_stack, evaluation_imaginary_term_count,
            integer_add(local_stack, division.quotient, &INT(1, -)));
    }
    else
    {
        evaluation_imaginary_term_count = integer_add(local_stack, division.quotient,
            integer_multiply(local_stack, output_stack, division.quotient, division.quotient));
        evaluation_real_term_count =
            integer_add(local_stack, evaluation_imaginary_term_count, division.quotient);
    }
    struct Rational*real_interval_scale = rational_magnitude(local_stack, a->coefficients[0]);
    struct Rational*imaginary_interval_scale = &rational_one;
    struct Rational*real_max_magnitude =
        number_estimate_max_magnitude(number_rational_real_part_estimate, pool_set, local_stack,
            output_stack, argument);
    struct Rational*imaginary_max_magnitude =
        number_estimate_max_magnitude(number_rational_imaginary_part_estimate, pool_set,
            local_stack, output_stack, argument);
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        struct Rational*coefficient_magnitude = rational_magnitude(local_stack, a->coefficients[i]);
        struct Rational*real_magnitude_power = rational_exponentiate(local_stack, output_stack,
            real_max_magnitude, integer_from_size_t(local_stack, i - 1));
        real_interval_scale = rational_max(output_stack, local_stack, real_interval_scale,
            rational_doubled(local_stack, output_stack,
                rational_multiply(local_stack, output_stack, coefficient_magnitude,
                    real_magnitude_power)));
        struct Rational*imaginary_magnitude_power = &rational_one;
        struct Integer*active_term_count = evaluation_imaginary_term_count;
        struct Integer*inactive_term_count = evaluation_real_term_count;
        for (size_t j = 1; j < i; ++j)
        {
            struct Rational*next_imaginary_magnitude_power = rational_multiply(local_stack,
                output_stack, imaginary_magnitude_power, imaginary_max_magnitude);
            struct Rational*shared_scale_component = rational_doubled(local_stack, output_stack,
                rational_integer_multiply(local_stack, output_stack,
                    rational_integer_multiply(local_stack, output_stack,
                        rational_multiply(local_stack, output_stack, coefficient_magnitude,
                            rational_add(local_stack, output_stack, &rational_one,
                                rational_add(local_stack, output_stack, real_magnitude_power,
                                    next_imaginary_magnitude_power))),
                        n_choose_k(local_stack, output_stack, integer_from_size_t(local_stack, i),
                            integer_from_size_t(local_stack, j))), active_term_count));
            real_magnitude_power = rational_divide(local_stack, output_stack, real_magnitude_power,
                real_max_magnitude);
            real_interval_scale = rational_max(output_stack, local_stack, real_interval_scale,
                rational_multiply(local_stack, output_stack, shared_scale_component,
                    real_magnitude_power));
            imaginary_interval_scale = rational_max(output_stack, local_stack,
                imaginary_interval_scale, rational_multiply(local_stack, output_stack,
                    shared_scale_component, imaginary_magnitude_power));
            imaginary_magnitude_power = next_imaginary_magnitude_power;
            POINTER_SWAP(active_term_count, inactive_term_count);
        }
        imaginary_interval_scale = rational_max(output_stack, local_stack, imaginary_interval_scale,
            rational_doubled(local_stack, output_stack,
                rational_multiply(local_stack, output_stack, imaginary_magnitude_power,
                    rational_integer_multiply(local_stack, output_stack, coefficient_magnitude,
                        active_term_count))));
    }
    interval_size = rational_place_value(local_stack, output_stack, interval_size);
    struct FloatInterval*argument_real_part_estimate =
        number_float_real_part_estimate(pool_set, local_stack, output_stack, argument,
            rational_divide(local_stack, output_stack, interval_size, real_interval_scale));
    struct FloatInterval*argument_imaginary_part_estimate =
        number_float_imaginary_part_estimate(pool_set, local_stack, output_stack, argument,
            rational_divide(local_stack, output_stack, interval_size, imaginary_interval_scale));
    struct RationalInterval evaluation_real_part = { a->coefficients[0], a->coefficients[0] };
    struct RationalInterval evaluation_imaginary_part = { &rational_zero, &rational_zero };
    struct Rational*real_part_min =
        float_to_rational(local_stack, output_stack, argument_real_part_estimate->min);
    struct Rational*real_part_max =
        float_to_rational(local_stack, output_stack, argument_real_part_estimate->max);
    struct Rational*imaginary_part_min =
        float_to_rational(local_stack, output_stack, argument_imaginary_part_estimate->min);
    struct Rational*imaginary_part_max =
        float_to_rational(local_stack, output_stack, argument_imaginary_part_estimate->max);
    for (size_t i = 1; i < a->coefficient_count; ++i)
    {
        struct Integer*term_degree = integer_from_size_t(local_stack, i);
        struct Rational*real_part_min_power =
            rational_exponentiate(local_stack, output_stack, real_part_min, term_degree);
        struct Rational*real_part_max_power =
            rational_exponentiate(local_stack, output_stack, real_part_max, term_degree);
        struct Rational*imaginary_part_min_power = &rational_one;
        struct Rational*imaginary_part_max_power = &rational_one;
        for (size_t j = 0; j <= i; ++j)
        {
            struct Rational*coefficient = rational_integer_multiply(local_stack, output_stack,
                a->coefficients[i], n_choose_k(local_stack, output_stack, term_degree,
                    integer_from_size_t(local_stack, j)));
            struct Rational*term_bound_candidate_a =
                rational_multiply(local_stack, output_stack, coefficient, real_part_min_power);
            struct Rational*term_bound_candidate_b =
                rational_multiply(local_stack, output_stack, coefficient, real_part_max_power);
            struct Rational*term_bound_candidate_c = rational_multiply(local_stack, output_stack,
                term_bound_candidate_a, imaginary_part_min_power);
            struct Rational*term_bound_candidate_d = rational_multiply(local_stack, output_stack,
                term_bound_candidate_b, imaginary_part_min_power);
            term_bound_candidate_a = rational_multiply(local_stack, output_stack,
                term_bound_candidate_a, imaginary_part_max_power);
            term_bound_candidate_b = rational_multiply(local_stack, output_stack,
                term_bound_candidate_b, imaginary_part_max_power);
            struct RationalInterval term_interval =
                { term_bound_candidate_a, term_bound_candidate_a };
            rational_interval_expand_bounds(output_stack, local_stack, &term_interval,
                term_bound_candidate_b);
            rational_interval_expand_bounds(output_stack, local_stack, &term_interval,
                term_bound_candidate_c);
            rational_interval_expand_bounds(output_stack, local_stack, &term_interval,
                term_bound_candidate_d);
            size_t degree_remainder = j % 4;
            switch (degree_remainder)
            {
            case 0:
                rational_add(local_stack, output_stack, evaluation_real_part.min,
                    term_interval.min);
                rational_add(local_stack, output_stack, evaluation_real_part.max,
                    term_interval.max);
                break;
            case 1:
                rational_add(local_stack, output_stack, evaluation_imaginary_part.min,
                    term_interval.min);
                rational_add(local_stack, output_stack, evaluation_imaginary_part.max,
                    term_interval.max);
                break;
            case 2:
                rational_subtract(local_stack, output_stack, evaluation_real_part.min,
                    term_interval.max);
                rational_subtract(local_stack, output_stack, evaluation_real_part.max,
                    term_interval.min);
                break;
            case 3:
                rational_subtract(local_stack, output_stack, evaluation_imaginary_part.min,
                    term_interval.max);
                rational_subtract(local_stack, output_stack, evaluation_imaginary_part.max,
                    term_interval.min);
                break;
            }
            real_part_min_power =
                rational_divide(local_stack, output_stack, real_part_min_power, real_part_min);
            real_part_max_power =
                rational_divide(local_stack, output_stack, real_part_max_power, real_part_max);
            imaginary_part_min_power = rational_multiply(local_stack, output_stack,
                imaginary_part_min_power, imaginary_part_min);
            imaginary_part_max_power = rational_multiply(local_stack, output_stack,
                imaginary_part_max_power, imaginary_part_max);
        }
    }
    rational_interval_to_float_interval(output_stack, local_stack, out->real_part_estimate,
        &evaluation_real_part, interval_size);
    rational_interval_to_float_interval(output_stack, local_stack, out->imaginary_part_estimate,
        &evaluation_imaginary_part, interval_size);
    local_stack->cursor = local_stack_savepoint;
}