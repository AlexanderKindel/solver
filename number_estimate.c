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
    struct Stack*output_stack, struct Stack*local_stack, struct FloatInterval*out,
    struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*rational_estimate_interval_size =
        rational_place_value(local_stack, output_stack, interval_size);
    struct RationalInterval rational_estimate;
    get_rational_estimate(local_stack, output_stack, &rational_estimate, a,
        rational_estimate_interval_size);
    rational_interval_to_float_interval(output_stack, local_stack, out, &rational_estimate,
        rational_estimate_interval_size);
    local_stack->cursor = local_stack_savepoint;
}

struct Rational*number_estimate_max_magnitude(rational_estimate_getter get_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval estimate;
    get_estimate(local_stack, output_stack, &estimate, a, &rational_one);
    struct Rational*out = rational_max(output_stack, local_stack,
        rational_magnitude(local_stack, estimate.min),
        rational_magnitude(local_stack, estimate.max));
    out = rational_copy(output_stack, out);
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_argument_cosine_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(local_stack, output_stack, &argument_estimate, a,
        argument_estimate_interval_size);
    struct RationalInterval argument_min_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_min_cosine, argument_estimate.min,
        argument_estimate_interval_size);
    struct RationalInterval argument_max_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_max_cosine, argument_estimate.max,
        argument_estimate_interval_size);
    out->max = rational_copy(output_stack,
        rational_max(output_stack, local_stack, argument_min_cosine.max, argument_max_cosine.max));
    pi_shrink_interval_to_one_side_of_value(argument_estimate.min);
    pi_shrink_interval_to_one_side_of_value(argument_estimate.max);
    struct Rational*pi_estimate_max =
        rational_add(local_stack, output_stack, pi_estimate_min, pi_interval_size);
    if (rational_compare(output_stack, local_stack, argument_estimate.min, pi_estimate_min) <= 0 &&
        rational_compare(output_stack, local_stack, pi_estimate_max, argument_estimate.max) <= 0)
    {
        out->min->numerator = integer_initialize(output_stack, 1, -1);
        out->min->denominator = &one;
    }
    else
    {
        out->min = rational_copy(output_stack, rational_max(output_stack, local_stack,
            argument_min_cosine.min, argument_max_cosine.min));
        out->min = rational_copy(output_stack, out->min);
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_argument_sine_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(local_stack, output_stack, &argument_estimate, a,
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
    pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.max);
    struct Rational*pi_estimate_max =
        rational_add(local_stack, output_stack, pi_estimate_min, pi_interval_size);
    if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
        pi_estimate_min) <= 0 && rational_compare(output_stack, local_stack, pi_estimate_max,
            argument_estimate_multiple.max) <= 0)
    {
        out->max = &rational_one;
        out->min = rational_copy(output_stack,
            rational_max(output_stack, local_stack, argument_min_sine.min, argument_max_sine.min));
    }
    else
    {
        argument_estimate_multiple.min = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.min, &INT(3, +));
        argument_estimate_multiple.max = rational_integer_divide(local_stack, output_stack,
            argument_estimate_multiple.max, &INT(3, +));
        pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.min);
        pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.max);
        pi_estimate_max =
            rational_add(local_stack, output_stack, pi_estimate_min, pi_interval_size);
        if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
            pi_estimate_min) <= 0 && rational_compare(output_stack, local_stack, pi_estimate_max,
                argument_estimate_multiple.max) <= 0)
        {
            out->min->numerator = integer_initialize(output_stack, 1, -1);
            out->min->denominator = &one;
            out->max = rational_copy(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
        else
        {
            out->min = rational_copy(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.min, argument_max_sine.min));
            out->max = rational_copy(output_stack, rational_max(output_stack, local_stack,
                argument_min_sine.max, argument_max_sine.max));
        }
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_rectangular_part_from_polar_form(void(trig_function)(struct Stack*, struct Stack*,
    struct RationalInterval*, struct Number*, struct Rational*),
    struct Stack*output_stack, struct Stack*local_stack, struct RationalInterval*out,
    struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval magnitude_estimate;
    number_rational_magnitude_estimate(local_stack, output_stack, &magnitude_estimate, a,
        &rational_one);
    struct Rational*factor_interval_size = rational_divide(local_stack, output_stack, interval_size,
        rational_integer_add(local_stack, output_stack, magnitude_estimate.max, &INT(2, +)));
    number_rational_magnitude_estimate(local_stack, output_stack, &magnitude_estimate, a,
        factor_interval_size);
    struct RationalInterval trig_value;
    trig_function(local_stack, output_stack, &trig_value, a, factor_interval_size);
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
    rational_estimate_getter get_factor_b_estimate, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*factor_a, struct Number*factor_b,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*out = rational_divide(output_stack, local_stack, interval_size,
        rational_add(local_stack, output_stack,
            number_estimate_max_magnitude(get_factor_a_estimate, local_stack, output_stack,
                factor_a),
            rational_integer_add(local_stack, output_stack,
                number_estimate_max_magnitude(get_factor_b_estimate, local_stack, output_stack,
                    factor_b), &one)));
    local_stack->cursor = local_stack_savepoint;
    return out;
}

void number_estimate_part_sum(float_estimate_getter get_term_estimate, struct Stack*output_stack,
    struct Stack*local_stack, struct FloatInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*term_estimate_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(2, +));
    struct FloatInterval left_term_estimate;
    get_term_estimate(local_stack, output_stack, &left_term_estimate, a->left,
        term_estimate_interval_size);
    struct FloatInterval right_term_estimate;
    get_term_estimate(local_stack, output_stack, &right_term_estimate, a->right,
        term_estimate_interval_size);
    out->min =
        float_add(output_stack, local_stack, left_term_estimate.min, right_term_estimate.min);
    out->max =
        float_add(output_stack, local_stack, left_term_estimate.max, right_term_estimate.max);
    local_stack->cursor = local_stack_savepoint;
}

void number_float_real_part_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
    {
        void*local_stack_savepoint = local_stack->cursor;
        number_float_magnitude_estimate(local_stack, output_stack, out, a, interval_size);
        if (a->value.numerator->sign < 0)
        {
            out->min = float_negative(output_stack, out->min);
            out->max = float_negative(output_stack, out->max);
        }
        else
        {
            out->min = float_copy(output_stack, out->min);
            out->max = float_copy(output_stack, out->max);
        }
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '^':
        number_float_estimate_from_rational(number_rational_real_part_estimate, output_stack,
            local_stack, out, a, interval_size);
        return;
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        interval_size = rational_half(local_stack, output_stack, interval_size);
        struct Rational*real_part_interval_size =
            factor_estimate_interval_size(number_rational_real_part_estimate,
                number_rational_real_part_estimate, local_stack, output_stack, a->left, a->right,
                interval_size);
        struct FloatInterval left_real_part;
        number_float_real_part_estimate(local_stack, output_stack, &left_real_part, a->left,
            real_part_interval_size);
        struct FloatInterval right_real_part;
        number_float_real_part_estimate(local_stack, output_stack, &right_real_part, a->left,
            real_part_interval_size);
        struct FloatInterval real_part_product;
        float_interval_multiply(local_stack, output_stack, &real_part_product, &left_real_part,
            &right_real_part);
        struct Rational*imaginary_part_interval_size =
            factor_estimate_interval_size(number_rational_imaginary_part_estimate,
                number_rational_imaginary_part_estimate, local_stack, output_stack, a->left,
                a->right, interval_size);
        struct FloatInterval left_imaginary_part;
        number_float_imaginary_part_estimate(local_stack, output_stack, &left_imaginary_part,
            a->left, imaginary_part_interval_size);
        struct FloatInterval right_imaginary_part;
        number_float_imaginary_part_estimate(local_stack, output_stack, &right_imaginary_part,
            a->right, imaginary_part_interval_size);
        struct FloatInterval imaginary_part_product;
        float_interval_multiply(local_stack, output_stack, &imaginary_part_product,
            &left_imaginary_part, &right_imaginary_part);
        out->min = float_subtract(output_stack, local_stack, real_part_product.min,
            imaginary_part_product.max);
        out->max = float_subtract(output_stack, local_stack, real_part_product.max,
            imaginary_part_product.min);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '+':
        number_estimate_part_sum(number_float_real_part_estimate, output_stack, local_stack, out, a,
            interval_size);
    }
}

void number_rational_real_part_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_cosine_estimate, output_stack,
            local_stack, out, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval float_estimate;
        number_float_real_part_estimate(local_stack, output_stack, &float_estimate, a,
            interval_size);
        out->min = float_to_rational(output_stack, local_stack, float_estimate.min);
        out->max = float_to_rational(output_stack, local_stack, float_estimate.max);
        local_stack->cursor = local_stack_savepoint;
    }
}

void number_float_imaginary_part_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        out->min = &float_zero;
        out->max = &float_zero;
        return;
    case '^':
        number_float_estimate_from_rational(number_rational_imaginary_part_estimate, output_stack,
            local_stack, out, a, interval_size);
        return;
    case '*':
    {
        void*local_stack_savepoint = local_stack->cursor;
        interval_size = rational_half(local_stack, output_stack, interval_size);
        struct Rational*term_a_interval_size =
            factor_estimate_interval_size(number_rational_real_part_estimate,
                number_rational_imaginary_part_estimate, local_stack, output_stack, a->left,
                a->right, interval_size);
        struct FloatInterval left_real_part;
        number_float_real_part_estimate(local_stack, output_stack, &left_real_part, a->left,
            term_a_interval_size);
        struct FloatInterval right_imaginary_part;
        number_float_imaginary_part_estimate(local_stack, output_stack, &right_imaginary_part,
            a->right, term_a_interval_size);
        struct FloatInterval term_a;
        float_interval_multiply(local_stack, output_stack, &term_a, &left_real_part,
            &right_imaginary_part);
        struct Rational*term_b_interval_size =
            factor_estimate_interval_size(number_rational_imaginary_part_estimate,
                number_rational_real_part_estimate, local_stack, output_stack, a->left, a->right,
                interval_size);
        struct FloatInterval left_imaginary_part;
        number_float_imaginary_part_estimate(local_stack, output_stack, &left_imaginary_part,
            a->left, term_b_interval_size);
        struct FloatInterval right_real_part;
        number_float_real_part_estimate(local_stack, output_stack, &right_real_part, a->right,
            term_b_interval_size);
        struct FloatInterval term_b;
        float_interval_multiply(local_stack, output_stack, &term_b, &left_imaginary_part,
            &right_real_part);
        out->min = float_add(output_stack, local_stack, term_a.min, term_b.min);
        out->max = float_add(output_stack, local_stack, term_a.max, term_b.max);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '+':
        number_estimate_part_sum(number_float_imaginary_part_estimate, output_stack, local_stack,
            out, a, interval_size);
    }
}

void number_rational_imaginary_part_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '^')
    {
        number_rectangular_part_from_polar_form(number_argument_sine_estimate, output_stack,
            local_stack, out, a, interval_size);
    }
    else
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct FloatInterval float_estimate;
        number_float_imaginary_part_estimate(local_stack, output_stack, &float_estimate, a,
            interval_size);
        out->min = float_to_rational(output_stack, local_stack, float_estimate.min);
        out->max = float_to_rational(output_stack, local_stack, float_estimate.max);
        local_stack->cursor = local_stack_savepoint;
    }
}

void number_rectangular_part_estimate_for_magnitude_estimate(
    float_estimate_getter get_part_float_estimate,
    rational_estimate_getter get_part_rational_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct FloatInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    get_part_float_estimate(output_stack, local_stack, out, a,
        rational_divide(local_stack, output_stack, interval_size,
            rational_integer_add(local_stack, output_stack,
                rational_integer_multiply(local_stack, output_stack,
                    number_estimate_max_magnitude(get_part_rational_estimate, local_stack,
                        output_stack, a), &INT(4, +)), &INT(2, +))));
    local_stack->cursor = local_stack_savepoint;
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

void number_float_magnitude_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    switch (a->operation)
    {
    case 'r':
    {
        rational_float_estimate(output_stack, local_stack, out, &a->value, interval_size);
        return;
    }
    case '^':
    {
        struct FloatInterval radicand_magnitude_estimate;
        number_float_magnitude_estimate(local_stack, output_stack, &radicand_magnitude_estimate,
            a->left, rational_integer_divide(local_stack, output_stack,
                rational_exponentiate(local_stack, output_stack, interval_size,
                    a->right->value.denominator), &INT(3, +)));
        struct Rational*bound_interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct Float*unused_estimate_bound;
        float_estimate_root(local_stack, output_stack, &out->min, &unused_estimate_bound,
            radicand_magnitude_estimate.min, bound_interval_size, a->right->value.denominator);
        float_estimate_root(local_stack, output_stack, &unused_estimate_bound, &out->max,
            radicand_magnitude_estimate.max, bound_interval_size, a->right->value.denominator);
        out->min = float_copy(output_stack, out->min);
        out->max = float_copy(output_stack, out->max);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '*':
    {
        struct RationalInterval left_factor_rational_magnitude_estimate;
        number_rational_magnitude_estimate(local_stack, output_stack,
            &left_factor_rational_magnitude_estimate, a->left, &rational_one);
        struct RationalInterval right_factor_rational_magnitude_estimate;
        number_rational_magnitude_estimate(local_stack, output_stack,
            &right_factor_rational_magnitude_estimate, a->right, &rational_one);
        interval_size = rational_divide(local_stack, output_stack, interval_size,
            rational_multiply(local_stack, output_stack,
                rational_add(local_stack, output_stack, left_factor_rational_magnitude_estimate.max,
                    &rational_one),
                rational_add(local_stack, output_stack,
                    right_factor_rational_magnitude_estimate.max, &rational_one)));
        struct FloatInterval left_factor_float_magnitude_estimate;
        number_float_magnitude_estimate(local_stack, output_stack,
            &left_factor_float_magnitude_estimate, a->left, interval_size);
        struct FloatInterval right_factor_float_magnitude_estimate;
        number_float_magnitude_estimate(local_stack, output_stack,
            &right_factor_float_magnitude_estimate, a->right, interval_size);
        out->min = float_multiply(output_stack, local_stack,
            left_factor_float_magnitude_estimate.min, right_factor_float_magnitude_estimate.min);
        out->max = float_multiply(output_stack, local_stack,
            left_factor_float_magnitude_estimate.max, right_factor_float_magnitude_estimate.max);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '+':
    {
        interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct FloatInterval real_part_estimate;
        number_rectangular_part_estimate_for_magnitude_estimate(number_float_real_part_estimate,
            number_rational_real_part_estimate, local_stack, output_stack, &real_part_estimate, a,
            interval_size);
        struct FloatInterval imaginary_part_estimate;
        number_rectangular_part_estimate_for_magnitude_estimate(
            number_float_imaginary_part_estimate, number_rational_imaginary_part_estimate,
            local_stack, output_stack, &imaginary_part_estimate, a, interval_size);
        struct Float*unused_estimate_bound;
        number_magnitude_estimate_bound_interval(local_stack, output_stack, &out->min,
            &unused_estimate_bound, real_part_estimate.min, imaginary_part_estimate.min,
            interval_size);
        number_magnitude_estimate_bound_interval(local_stack, output_stack, &unused_estimate_bound,
            &out->max, real_part_estimate.max, imaginary_part_estimate.max, interval_size);
        out->min = float_copy(output_stack, out->min);
        out->max = float_copy(output_stack, out->max);
        local_stack->cursor = local_stack_savepoint;
    }
    }
}

void number_rational_magnitude_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct FloatInterval float_estimate;
    number_float_magnitude_estimate(local_stack, output_stack, &float_estimate, a, interval_size);
    float_interval_to_rational_interval(output_stack, local_stack, out, &float_estimate);
    local_stack->cursor = local_stack_savepoint;
}

void float_argument_estimate_for_part_sum(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    number_float_argument_estimate(output_stack, local_stack, out, a, interval_size);
}

void number_float_argument_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct FloatInterval*out, struct Number*a, struct Rational*interval_size)
{
    if (a->operation == '*')
    {
        number_estimate_part_sum(float_argument_estimate_for_part_sum, output_stack, local_stack,
            out, a, interval_size);
    }
    else
    {
        number_float_estimate_from_rational(number_rational_argument_estimate, output_stack,
            local_stack, out, a, interval_size);
    }
}

void number_rational_argument_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    switch (a->operation)
    {
    case 'r':
        if (a->value.numerator->sign < 0)
        {
            pi_estimate(interval_size);
            out->min = pi_estimate_min;
            out->max = rational_add(output_stack, local_stack, pi_estimate_min, pi_interval_size);
        }
        else
        {
            out->min = &rational_zero;
            out->max = &rational_zero;
        }
        return;
    case '^':
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval radicand_rational_argument_estimate;
        number_rational_argument_estimate(local_stack, output_stack,
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
        struct FloatInterval float_estimate;
        number_float_argument_estimate(local_stack, output_stack, &float_estimate, a,
            interval_size);
        float_interval_to_rational_interval(output_stack, local_stack, out, &float_estimate);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    case '+':
        crash("number_rational_argument_estimate case not yet implemented.");
    }
}

void number_rectangular_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RectangularEstimate*out, struct Number*a, struct Rational*interval_size)
{
    number_float_real_part_estimate(output_stack, local_stack, out->real_part_estimate, a,
        interval_size);
    number_float_imaginary_part_estimate(output_stack, local_stack, out->imaginary_part_estimate, a,
        interval_size);
}

bool rectangular_estimates_are_disjoint(struct Stack*stack_a, struct Stack*stack_b,
    struct RectangularEstimate*a, struct RectangularEstimate*b)
{
    return float_intervals_are_disjoint(stack_a, stack_b, a->real_part_estimate,
            b->real_part_estimate) ||
        float_intervals_are_disjoint(stack_a, stack_b, a->imaginary_part_estimate,
            b->imaginary_part_estimate);
}

void rational_polynomial_estimate_evaluation(struct Stack*output_stack, struct Stack*local_stack,
    struct RectangularEstimate*out, struct RationalPolynomial*a, struct Number*argument,
    struct Rational*interval_size)
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
        rational_float_estimate(output_stack, local_stack, out->real_part_estimate,
            a->coefficients[0], interval_size);
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
        number_estimate_max_magnitude(number_rational_real_part_estimate, local_stack, output_stack,
            argument);
    struct Rational*imaginary_max_magnitude =
        number_estimate_max_magnitude(number_rational_imaginary_part_estimate, local_stack,
            output_stack, argument);
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
    struct FloatInterval argument_real_part_estimate;
    number_float_real_part_estimate(local_stack, output_stack, &argument_real_part_estimate,
        argument, rational_divide(local_stack, output_stack, interval_size, real_interval_scale));
    struct FloatInterval argument_imaginary_part_estimate;
    number_float_imaginary_part_estimate(local_stack, output_stack,
        &argument_imaginary_part_estimate, argument,
        rational_divide(local_stack, output_stack, interval_size, imaginary_interval_scale));
    struct RationalInterval evaluation_real_part = { a->coefficients[0], a->coefficients[0] };
    struct RationalInterval evaluation_imaginary_part = { &rational_zero, &rational_zero };
    struct Rational*real_part_min =
        float_to_rational(local_stack, output_stack, argument_real_part_estimate.min);
    struct Rational*real_part_max =
        float_to_rational(local_stack, output_stack, argument_real_part_estimate.max);
    struct Rational*imaginary_part_min =
        float_to_rational(local_stack, output_stack, argument_imaginary_part_estimate.min);
    struct Rational*imaginary_part_max =
        float_to_rational(local_stack, output_stack, argument_imaginary_part_estimate.max);
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