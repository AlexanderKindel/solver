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

void number_refine_rational_estimate_interval(rational_estimate_getter get_rational_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct RationalInterval*out,
    struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    if (rational_compare(output_stack, local_stack, interval_size,
        rational_subtract(local_stack, output_stack, out->max, out->min)) < 0)
    {
        get_rational_estimate(output_stack, local_stack, out, a, interval_size);
    }
    local_stack->cursor = local_stack_savepoint;
}

void number_halve_rational_estimate_interval(rational_estimate_getter get_rational_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct RationalInterval*out,
    struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    get_rational_estimate(output_stack, local_stack, out, a,
        rational_half(local_stack, output_stack,
            rational_subtract(local_stack, output_stack, out->max, out->min)));
    local_stack->cursor = local_stack_savepoint;
}

void number_halve_float_estimate_interval(float_estimate_getter get_float_estimate,
    struct Stack*output_stack, struct Stack*local_stack, struct FloatInterval*out, struct Number*a)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Float*current_interval_size =
        float_subtract(local_stack, output_stack, out->max, out->min);
    get_float_estimate(output_stack, local_stack, out, a,
        float_to_rational(local_stack, output_stack,
            &(struct Float){current_interval_size->significand,
                integer_add(local_stack, current_interval_size->significand, &INT(1, -))}));
    local_stack->cursor = local_stack_savepoint;
}

void number_rectangular_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RectangularEstimate*out, struct Number*a, struct Rational*interval_size)
{
    number_float_real_part_estimate(output_stack, local_stack, out->real_part_estimate, a,
        interval_size);
    number_float_imaginary_part_estimate(output_stack, local_stack, out->imaginary_part_estimate, a,
        interval_size);
}

void number_rectangular_estimate_halve_dimensions(struct Stack*output_stack,
    struct Stack*local_stack, struct RectangularEstimate*out, struct Number*a)
{
    number_halve_float_estimate_interval(number_float_real_part_estimate, output_stack, local_stack,
        out->real_part_estimate, a);
    number_halve_float_estimate_interval(number_float_imaginary_part_estimate, output_stack,
        local_stack, out->imaginary_part_estimate, a);
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

void number_argument_cosine_estimate(struct Stack*output_stack, struct Stack*local_stack,
    struct RationalInterval*out, struct Number*a, struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct Rational*argument_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(local_stack, output_stack, &argument_estimate, a,
        argument_interval_size);
    struct RationalInterval argument_min_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_min_cosine, argument_estimate.min,
        argument_interval_size);
    struct RationalInterval argument_max_cosine;
    rational_estimate_cosine(local_stack, output_stack, &argument_max_cosine, argument_estimate.max,
        argument_interval_size);
    out->max = rational_copy(output_stack, rational_max(output_stack, local_stack,
        argument_min_cosine.max, argument_max_cosine.max));
    pi_shrink_interval_to_one_side_of_value(argument_estimate.min);
    pi_shrink_interval_to_one_side_of_value(argument_estimate.max);
    if (rational_compare(output_stack, local_stack, argument_estimate.min, pi.min) <= 0 &&
        rational_compare(output_stack, local_stack, pi.max, argument_estimate.max) <= 0)
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
    struct Rational*argument_interval_size =
        rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
    struct RationalInterval argument_estimate;
    number_rational_argument_estimate(local_stack, output_stack, &argument_estimate, a,
        argument_interval_size);
    struct RationalInterval argument_min_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_min_sine, argument_estimate.min,
        argument_interval_size);
    struct RationalInterval argument_max_sine;
    rational_estimate_sine(local_stack, output_stack, &argument_max_sine, argument_estimate.max,
        argument_interval_size);
    struct RationalInterval argument_estimate_multiple =
    { rational_doubled(local_stack, output_stack, argument_estimate.min),
    rational_doubled(local_stack, output_stack, argument_estimate.max) };
    pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.min);
    pi_shrink_interval_to_one_side_of_value(argument_estimate_multiple.max);
    if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min, pi.min) <= 0 &&
        rational_compare(output_stack, local_stack, pi.max, argument_estimate_multiple.max) <= 0)
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
        if (rational_compare(output_stack, local_stack, argument_estimate_multiple.min,
            pi.min) <= 0 && rational_compare(output_stack, local_stack, pi.max,
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
    number_refine_rational_estimate_interval(number_rational_magnitude_estimate, local_stack,
        output_stack, &magnitude_estimate, a, factor_interval_size);
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

struct FloatInterval*number_estimate_product(struct EstimateGetters*factor_a_estimate_getters,
    struct EstimateGetters*factor_b_estimate_getters, struct Stack*output_stack,
    struct Stack*local_stack, struct Number*factor_a, struct Number*factor_b,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval rational_factor_a_estimate;
    factor_a_estimate_getters->rational(local_stack, output_stack, &rational_factor_a_estimate,
        factor_a, interval_size);
    struct RationalInterval rational_factor_b_estimate;
    factor_b_estimate_getters->rational(local_stack, output_stack, &rational_factor_b_estimate,
        factor_b, interval_size);
    interval_size = rational_interval_factor_interval_size(local_stack, output_stack,
        &rational_factor_a_estimate, &rational_factor_b_estimate, interval_size);
    struct FloatInterval factor_a_estimate;
    factor_a_estimate_getters->fl(output_stack, local_stack, &factor_a_estimate, factor_a,
        interval_size);
    struct FloatInterval factor_b_estimate;
    factor_b_estimate_getters->fl(output_stack, local_stack, &factor_b_estimate, factor_b,
        interval_size);
    struct FloatInterval*out =
        float_interval_multiply(output_stack, local_stack, &factor_a_estimate, &factor_b_estimate);
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
        struct FloatInterval*real_part_product = number_estimate_product(&real_estimate_getters,
            &real_estimate_getters, local_stack, output_stack, a->left, a->right, interval_size);
        struct FloatInterval*imaginary_part_product =
            number_estimate_product(&imaginary_estimate_getters, &imaginary_estimate_getters,
                local_stack, output_stack, a->left, a->right, interval_size);
        out->min = float_subtract(output_stack, local_stack, real_part_product->min,
            imaginary_part_product->max);
        out->max = float_subtract(output_stack, local_stack, real_part_product->max,
            imaginary_part_product->min);
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
        struct FloatInterval*term_a =
            number_estimate_product(&real_estimate_getters, &imaginary_estimate_getters,
                local_stack, output_stack, a->left, a->right, interval_size);
        struct FloatInterval*term_b = number_estimate_product(&imaginary_estimate_getters,
            &real_estimate_getters, local_stack, output_stack, a->left, a->right, interval_size);
        out->min = float_add(output_stack, local_stack, term_a->min, term_b->min);
        out->max = float_add(output_stack, local_stack, term_a->max, term_b->max);
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
    struct EstimateGetters*part_estimate_getters, struct Stack*output_stack,
    struct Stack*local_stack, struct FloatInterval*out, struct Number*a,
    struct Rational*interval_size)
{
    void*local_stack_savepoint = local_stack->cursor;
    struct RationalInterval estimate;
    part_estimate_getters->rational(local_stack, output_stack, &estimate, a, &rational_one);
    part_estimate_getters->fl(output_stack, local_stack, out, a,
        rational_divide(local_stack, output_stack, interval_size,
            rational_integer_add(local_stack, output_stack,
                rational_integer_multiply(local_stack, output_stack,
                    rational_interval_max_magnitude(local_stack, output_stack, &estimate),
                    &INT(4, +)), &INT(2, +))));
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
        struct FloatInterval*product = number_estimate_product(&magnitude_estimate_getters,
            &magnitude_estimate_getters, local_stack, output_stack, a->left, a->right,
            interval_size);
        out->min = float_copy(output_stack, product->min);
        out->max = float_copy(output_stack, product->max);
        return;
    }
    case '+':
    {
        interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct FloatInterval real_part_estimate;
        number_rectangular_part_estimate_for_magnitude_estimate(&real_estimate_getters, local_stack,
            output_stack, &real_part_estimate, a, interval_size);
        struct FloatInterval imaginary_part_estimate;
        number_rectangular_part_estimate_for_magnitude_estimate(&imaginary_estimate_getters,
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
            out->min = pi.min;
            out->max = pi.max;
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
    {
        void*local_stack_savepoint = local_stack->cursor;
        struct RationalInterval real_part_estimate;
        number_rational_real_part_estimate(local_stack, output_stack, &real_part_estimate, a,
            &rational_one);
        struct RationalInterval imaginary_part_estimate;
        number_rational_imaginary_part_estimate(local_stack, output_stack,
            &imaginary_part_estimate, a, &rational_one);
        if (imaginary_part_estimate.min->numerator->sign !=
            imaginary_part_estimate.max->numerator->sign)
        {
            while (true)
            {
                struct Rational*imaginary_max_magnitude =
                    rational_interval_max_magnitude(local_stack, output_stack,
                        &imaginary_part_estimate);
                if (rational_polynomial_root_count_in_rectangle(local_stack, output_stack,
                    a->minimal_polynomial, &real_part_estimate,
                    &(struct RationalInterval){rational_negative(local_stack,
                        imaginary_max_magnitude), imaginary_max_magnitude}) == 1)
                {
                    while (real_part_estimate.min->numerator->sign ==
                        -real_part_estimate.max->numerator->sign)
                    {
                        number_halve_rational_estimate_interval(number_rational_real_part_estimate,
                            local_stack, output_stack, &real_part_estimate, a);
                    }
                    if (real_part_estimate.max->numerator->sign > 0)
                    {
                        out->min = &rational_zero;
                        out->max = &rational_zero;
                    }
                    else
                    {
                        pi_estimate(interval_size);
                        out->min = pi.min;
                        out->max = pi.max;
                    }
                    local_stack->cursor = local_stack_savepoint;
                    return;
                }
                number_halve_rational_estimate_interval(number_rational_imaginary_part_estimate,
                    local_stack, output_stack, &imaginary_part_estimate, a);
                if (imaginary_part_estimate.min->numerator->sign ==
                    imaginary_part_estimate.max->numerator->sign)
                {
                    break;
                }
                number_halve_rational_estimate_interval(number_rational_real_part_estimate,
                    local_stack, output_stack, &real_part_estimate, a);
            }
        }
        interval_size =
            rational_integer_divide(local_stack, output_stack, interval_size, &INT(3, +));
        struct Rational*imaginary_part_interval_size =
            rational_interval_factor_interval_size(local_stack, output_stack,
                &imaginary_part_estimate, &(struct RationalInterval){
                    rational_reciprocal(local_stack, real_part_estimate.max),
                    rational_reciprocal(local_stack, real_part_estimate.min)}, interval_size);
        number_refine_rational_estimate_interval(number_rational_imaginary_part_estimate,
            local_stack, output_stack, &imaginary_part_estimate, a, imaginary_part_interval_size);
        struct Rational*real_part_interval_size = rational_multiply(local_stack, output_stack,
            real_part_estimate.min, imaginary_part_interval_size);
        real_part_interval_size = rational_divide(local_stack, output_stack,
            rational_multiply(local_stack, output_stack, real_part_estimate.min,
                real_part_interval_size),
            rational_subtract(local_stack, output_stack, &rational_one, real_part_interval_size));
        number_refine_rational_estimate_interval(number_rational_real_part_estimate,
            local_stack, output_stack, &real_part_estimate, a, real_part_interval_size);
        rational_interval_estimate_atan2(output_stack, local_stack, out, &imaginary_part_estimate,
            &real_part_estimate, interval_size);
        local_stack->cursor = local_stack_savepoint;
        return;
    }
    }
}